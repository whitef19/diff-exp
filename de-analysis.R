#--------------------------------
# Differential expression analysis with SVA and limma-voom  
# by Frederique White
# 2022-11-24 
# 
# Usage : Rscript de-analysis.R -c path_to_config_file -o output_directory
#
# 
# Libraries
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(reshape2))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(require(ggrepel))
suppressMessages(library(ggpubr))

#--------------------------------

if(!exists("usage", mode="function")) source("de-functions.R")

prepare_config = function(config_path) 
{
	if (!(file.exists(config_path))) {
		usage("Error : config file not found")
	}
	
	config = read.csv(config_path, sep="\t", header=T, row.names=1)
	if(ncol(config) != 2) {
		stop("config file must have 3 columns: variable names, types and values")
	} 
	write.table(cbind(variable_name=rownames(config), config), file=paste0(out_path, "/config_file_used"), sep="\t", quote=F, row.names=F)
	colnames(config) <- c("variable_type", "variable_value")
	paths = data.frame(t(subset(config, variable_type == "file", "variable_value")))
	parameters = data.frame(t(subset(config, variable_type == "parameter", "variable_value")))

	# verify paths
	for (file in paths) {
		if (!is.na(file)){
			if (!(file.exists(file))) {
				stop(paste0("Error : ",file," file not found"))}
		}
	}

	return(list(files=paths, parameters=parameters))	
}

prepare_count_matrix = function(file_name) 
{
	## look if you have to skipped header rows
	skipped = 0
	res <- try(read.csv(file_name, sep="\t", check.names=F, row.names=1, skip=skipped))
	while (inherits(res, "try-error")) {
		skipped = skipped + 1
		res <- try(read.csv(file_name, sep="\t", check.names=F, row.names=1, skip=skipped))
	}
	
	count_df = read.csv(file_name, sep="\t", check.names=F, row.names=1, skip=skipped)
	sample_id = colnames(count_df)
	cat(paste0("=> initial number of genes: ", nrow(count_df), "\n"), file=log)
	
	## remove first column that contains description
	if (!is.numeric(count_df[1,1])){
		count_df[, 1] <- NULL
	}
	count_df <- data.frame(count_df)
	colnames(count_df) <- sample_id

	return(count_df)
}

prepare_design <- function(design_file, sample_file, params)
{
	design = read.csv(design_file, sep="\t", check.names=F, row.names=1)
	included_samples = read.csv(sample_file, header=F)$V1
	design = design[included_samples, ]

	full_model = paste0(" ~ ", params$full_model)
	null_model = paste0(" ~ ", params$null_model)

	#### put covariates as factors
	#included_covariates <- unlist(strsplit(params$model, "+", fixed=T))
	#for (cov in included_covariates) {
	#	if (cov != params$contrast) {
	#		if ( nrow(data.frame(table(design[, cov]))) < 5) {
	#			design[, cov] = factor(design[, cov])
	#		}
	#	}
	#}


	model_matrix = model.matrix(eval(parse(text=full_model)), data=design) 
	included_samples = intersect(included_samples, rownames(model_matrix))
	cat(paste0("=> number of included samples: ", length(included_samples),"\n"),file=log)

	design = design[included_samples, ]
	return(list(design=design, full_model=full_model, null_model=null_model, samples=included_samples))
}


gene_filtering = function(counts, norm_counts, design, params, map_scores_file)
{
	
	### Apply count thresholds
	cat(paste0("Gene filter thresholds:\ncounts: ", params$count_threshold,"\nsample fraction: ", params$sample_frac_threshold,"\n"), file=log)
	passed_genes <- apply(counts, 1, function(x) (sum(x >= as.numeric(params$count_threshold))/length(x)) >= as.numeric(params$sample_frac_threshold))
	
	if (exists("norm_counts", mode="list")) {
		cat(paste0("normalized counts: ", params$normalized_count_threshold,"\n"), file=log)
		passed_norm_thres <- apply(norm_counts, 1, function(x) (sum(x >= as.numeric(params$normalized_count_threshold))/length(x)) >= as.numeric(params$sample_frac_threshold))
		passed_genes <- (passed_genes & passed_norm_thres)
	}

	included_genes <- rownames(counts)[passed_genes]
	cat(paste0("=> passed genes: ", length(included_genes),"\n"), file=log)

	### Apply mappability threshold
	if (!is.na(map_scores_file)) {
		
		mappability_scores <- read.table(map_scores_file, header=F, col.names=c("gene_id","score"))
		high_mappability_genes <- subset(mappability_scores, score >= params$mappability_threshold)$gene_id
		included_genes <- intersect(high_mappability_genes, included_genes)
		cat(paste0("# mappability score threshold: ", params$mappability_threshold,"\n=> high mappability genes: ", length(included_genes),"\n"), file=log)
	}

	return(included_genes)
}


run_sva = function(counts, design, full_model, null_model){
	
	suppressMessages(library(edgeR))
	suppressMessages(library(limma))
	suppressMessages(library(isva))
	suppressMessages(library(SmartSVA))
	suppressMessages(library(corrplot))

	full_model_matrix = model.matrix(eval(parse(text=full_model)), data=design) 
	null_model_matrix = model.matrix(eval(parse(text=null_model)), data=design)
	
	voom = edger_norm(counts, full_model_matrix)$voom$E
	write.table(cbind("ID"=rownames(voom), voom), file=paste0(out_path, "/voom_counts.tsv" ), sep="\t", quote=F, row.names=F)

	numSVs = how_much_svs(voom, design, full_model, full_model_matrix, null_model_matrix)
	print(numSVs)
	print(num.sv(voom, null_model_matrix, method="leek"))

	#pValues = f.pvalue(voom,full_model_matrix, null_model_matrix)
	#qValues = p.adjust(pValues, method="BH")
	
	### run sva
	cat(paste0("# Compute surrugate variables\n", " Full model: ", full_model ,"\n", " Null model: ", null_model ,"\n"), file=log)
	SVObject <- SmartSVA::smartsva.cpp(voom, full_model_matrix, mod0=null_model_matrix, n.sv=numSVs, alpha=1, B=200, VERBOSE=F)

	### concatenate SVs to covariate dataframe 
	SVs <- as.data.frame(SVObject$sv)
	names(SVs) <- paste0("SV",1:ncol(SVs))
	design <- cbind(design, SVs)
	
	#modSv = cbind(full_model_matrix, SVObject$sv)
	#mod0Sv = cbind(null_model_matrix, SVObject$sv)
	#pValuesSv = f.pvalue(voom,modSv,mod0Sv)
	#qValuesSv = p.adjust(pValuesSv, method="BH")

	#df <- data.frame(ID=names(qValues), qvalue=qValues, qvaluesv=qValuesSv)
	#print(ncol(subset(df, qvalue<0.05)))
	#print(ncol(subset(df, qvaluesv<0.05)))

	correlation_matrix(design, numSVs)
	
	write.table(cbind("ID"=rownames(design), design), file=paste0(out_path,"/design_svs.tsv"), sep="\t", row.names=F, quote=F)
	
	residuals <- t(resid(lm(as.formula(paste0('t(voom) ', "~", paste0("SV",1:numSVs, collapse="+"))), data=design)))
	write.table(cbind("ID"=rownames(residuals), residuals), file=paste0(out_path,"/sv_residuals.tsv"), sep="\t", row.names=F, quote=F)

	full_model = paste0(full_model, "+", paste0("SV",1:numSVs, collapse="+"))

	return(list(design=design, model=full_model))
}




run_limma = function(counts, design, annotations, full_model)
{
	### create covariate matrix
	model_matrix = model.matrix(eval(parse(text=full_model)), data=design)

	### rerun voom with full model
	norm = edger_norm(counts, model_matrix)
	voom = norm$voom
	voom_counts = voom$E
	write.table(cbind("ID"=rownames(voom_counts), voom_counts), file=paste0(out_path, "/voom_counts.tsv" ), sep="\t", quote=F, row.names=F)

	full_model_but_contrast = paste0("~", paste0(unlist(strsplit(full_model, "+", fixed=T))[-1], collapse="+"))
	residuals <- t(resid(lm(as.formula(paste0('t(voom_counts) ', full_model_but_contrast)), data=design)))
	write.table(cbind("ID"=rownames(residuals), residuals), file=paste0(out_path,"/full_model_residuals.tsv"), sep="\t", row.names=F, quote=F)

	### run the regressions
	fit <- limma::lmFit(voom, model_matrix)
	contrast = colnames(fit$coefficients)[2]

	### run eBayes() to borrow information across genes
	ebayes_fit <- limma::eBayes(fit)
	save(ebayes_fit, file=paste0(out_path,"/ebayes_fit.RData"))

	### Extract limma results
	results <- limma::topTable(ebayes_fit, coef=contrast, number=Inf, sort.by="p")
	results <- cbind(ID=rownames(results), results)
	
	### Annotate limma results
	average_tmm <- data.frame(ID=rownames(norm$tmm), avgTMM=rowMeans(norm$tmm))
	results <- merge(results, average_tmm, by="ID", all.x=T)
	
	if (file.exists(annotations)) {
		
		gene_description = read.csv(annotations, sep="\t")
		results <- merge(gene_description, results, by="ID", all.y=T)
	}

	results <- results[order(results$P.Value), ]
	write.table(results, file=paste0(out_path,"/limma_results.tsv"), sep="\t", quote=F, row.names=F)
	
	return(results)
}

run_deseq <- function(counts, design, config, normalized_count_included, full_model)
{

	ddsFullCountTable = DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=design, design=as.formula(full_model))
	dds <- DESeq2::DESeq(ddsFullCountTable)
	print(DESeq2::resultsNames(dds))
	contrast = DESeq2::resultsNames(dds)[2]
	print(paste0("DESeq2 is using : ",contrast, " as contrast !"))
	deseq_counts <- DESeq2::counts(dds, normalized=TRUE)
	deseq_counts <- cbind(ID=rownames(deseq_counts), deseq_counts)
	write.table(deseq_counts, file=paste0(out_path, "/deseq_normalized_counts.tsv" ), sep="\t", quote=F, row.names=F)
	
	results <- data.frame(DESeq2::results(dds, name=contrast))
	results <- cbind(ID=rownames(results), results)	

	if (! is.na(config$paths$gene_annotations)) {
		gene_description = read.csv(config$paths$gene_annotations, sep="\t")
		results <- merge(gene_description, results, by="ID", all.y=T)
		colnames(results) <- c("ID",colnames(gene_description)[-1],"baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
	} else {
		colnames(results) <- c("ID","baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
	}

	results <- results[order(results$P.Value),]
	write.table(results, file=paste0(out_path,"/deseq_results.tsv"), sep="\t", quote=F, row.names=F)
	return(results)
}

create_figures = function(results, model, analysis, design)
{

	nbsvs = colnames(design)[ncol(design)]
	nbsamples = nrow(design)
	results$absFC = abs(results$logFC)
	logFC_threshold <- (mean(results$absFC) + 5 * sd(results$absFC))

	logFC_threshold = 0
	
	#top_hits <- subset(results, absFC > logFC_threshold)
	top_hits <- subset(results, P.Value < 0.05 )
	#results$color <- ifelse(results$absFC > logFC_threshold, "#22908C", "#F2BB05")
	#results$color <- ifelse(results$P.Value > 0.05, "black", results$color)
	results$color <- ifelse(results$logFC > 0, "green", "red")
	results$color <- ifelse(results$P.Value > 0.05, "grey", results$color)


	cat(paste0("logFC threshold mean + 5 SDs: ",logFC_threshold,  "\n","==> number of top hits: ",nrow(top_hits),"\n"), file=log)
	write.table(top_hits, file=paste0(out_path,"/top_", analysis, "_results.tsv"), sep="\t", quote=F, row.names=F)

	pdf(paste0(out_path, "/figures.",analysis,".pdf"), width=10, height=6)

	volcano <- ggplot(results, aes(x=logFC, y=-log10(P.Value), color=color)) + geom_point(alpha=0.5) + 
		ylab("-log10 nominal p-value") + xlab("log2 FC") + 
		xlim(-max(results$absFC), max(results$absFC)) +
		ggtitle(paste0("~ ", model, "\n + ", nbsvs,"\n n = ", nbsamples, " samples")) + 
		#geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), color="grey", linetype="dashed") +
		#geom_hline(yintercept=c(-log10(0.05), -log10(0.001)), color="grey",linetype="dashed") +
		geom_hline(yintercept=-log10(0.05), color="grey", linetype="dashed") +
		theme(plot.title=element_text(size=5), legend.position="none") + scale_color_identity()

	if (nrow(subset(results, absFC > logFC_threshold & P.Value < 0.05)) > 0){
		if ( "Name" %in% names(results)) {
			print(volcano + geom_text_repel(data =subset(results, absFC > logFC_threshold & P.Value < 0.05), aes(label=Name, size=1), show.legend=F, colour = "black"))
		} else {
			print(volcano + geom_text_repel(data=subset(results, absFC > logFC_threshold & P.Value < 0.05), aes(label=ID, size=1), show.legend=F, colour = "black"))
		}
	} else {print(volcano)}

	dev.off()

}

boxplots <- function(gene_name, logcpm){
	gene_name <- gsub("-", ".", gene_name)
	p1 <- ggplot(logcpm, aes(x=as.factor(contrast), y=as.numeric(logcpm[, gene_name]), group=as.factor(contrast), fill=as.factor(contrast))) + xlab("") + ylab(paste0("log2 CPM ", gene_name )) + geom_boxplot() + stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	print(p1)
}

scatterplot <- function(gene_name, logcpm){
	gene_name <- gsub("-", ".", gene_name)
	p1 <- ggplot(logcpm, aes(x=contrast, y=as.numeric(logcpm[, gene_name]))) + ylab(paste0("log2 CPM ", gene_name )) + geom_point() 
	print(p1)
}

#------
#-----------
#---------------------------
#-------------------
#----------------------------------

        #
        ###
        ####
##############     ======
###############     MAIN 
##############     ======
        ####
        ###
        #

#----------------------------------
#-------------------
#---------------------------
#-----------
#------
#           

ARG = commandArgs(trailingOnly = T)
if (length(ARG) < 2) { usage("missing arguments") }

## default arg values
out_path = "./"

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-c") {
		config_path = ARG[i+1]
	} else if (ARG[i] == "-o") {
		out_path = ARG[i+1]
	} else if (ARG[i] == "-h") {
		usage("")
	}
}

if (! dir.exists(out_path)) { dir.create(out_path) }
log = file(paste0(out_path,"/errors.o"), open="a")

## read config file and verify if all provided files exist
config = prepare_config(config_path)
parameters = config$parameters
files = config$files


## open count file
counts = prepare_count_matrix(files$read_count_matrix)
if (is.na(files$normalized_count_matrix)) { 
	normalized_counts = NA
} else {
	normalized_counts = prepare_count_matrix(files$normalized_count_matrix)
}


## read design file and list of sample to included in analysis
des = prepare_design(files$design, files$samples, parameters)
design = des$design
full_model = des$full_model
null_model = des$null_model
included_samples = des$samples


### filter count matrix
included_genes = gene_filtering(counts, normalized_counts, design, parameters, files$mappability_scores)
clean_count = counts[included_genes, included_samples]
write.table(cbind("ID"=rownames(clean_count), clean_count), file=paste0(out_path,"/clean_counts.tsv"), sep="\t", quote=F, row.names=F)


### run surrogate variable analysis
if (parameters$run_SVA) {

	sva = run_sva(clean_count, design, full_model, null_model)
	design = sva$design
	full_model = sva$model

} else if (parameters$number_surrogate_variables > 0) {

	full_model = paste0("~", full_model, "+", paste0("SV",1:parameters$number_surrogate_variables, collapse="+"))

} else {

	full_model = paste0("~", full_model)
}

cat(paste0("# Full Model: ", full_model ,"\n"), file=log)

### Differential expression analysis 

if (parameters$run_limma){

	results = run_limma(clean_count, design, files$gene_annotations, full_model)
	
	create_figures(results, full_model, "limma", design)

	#logcpm <- log2(cpm + 1)	
	#logcpm <- data.frame( cbind(ID=colnames(logcpm), design[colnames(logcpm), ], data.frame(t(logcpm)) ) )
	#logcpm$contrast = design[, parameters$contrast]
	#is_dichotomic <- nrow(data.frame(table(design[, parameters$contrast]))) == 2

	#pdf(paste0(out_path, "/top_genes.pdf"))
	#for(i in results[1:10, "ID"]) {
		#if (is_dichotomic) {
			#boxplots(i, logcpm)
		#} else {
			#scatterplot(i, logcpm)
		#}
	#}
	#dev.off()
}
if (parameters$run_deseq){
	print("Run DE analysis with DESeq2")
	results_deseq = run_deseq(count_df, design, config, normalized_count_included, full_model)
	create_figures(results_deseq, normalized_count_included, nb_svs,"deseq")

}












#violin_plot <- function(gene_name, gene_levels){
#	suppressMessages(library(ggpubr))
#	p1 <- ggplot(gene_levels, aes(x=categories, y=normalized, fill=categories)) + xlab("") + ylab(paste0("TPM ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
#	p2 <- ggplot(gene_levels, aes(x=categories, y=voom, fill=categories)) + xlab("") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
#	p3 <- ggplot(gene_levels, aes(x=categories, y=residuals, fill=categories)) + xlab("") + ylab(paste0("Full model residuals ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
#	print(ggarrange(p1, p2, p3, ncol=3))
#}
#
#scatter_plot <- function(gene_name, gene_levels){
#	suppressMessages(library(ggpubr))
#	p1 <- ggplot(gene_levels, aes(x=BMI, y=TPM)) + xlab("log2 BMI") + ylab(paste0("TPM ", gene_name )) + geom_point(alpha=0.5,color="#22908C") 
#	p2 <- ggplot(gene_levels, aes(x=BMI, y=voom)) + xlab("log2 BMI") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
#	p3 <- ggplot(gene_levels, aes(x=BMI, y=residuals)) + xlab("log2 BMI") + ylab(paste0("Full model residuals ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
#	print(ggarrange(p1, p2, p3, ncol=3))
#}

#contrast_is_dicho <- nrow(as.data.frame(table(covariates[, contrast]))) == 2
#top_genes <- results[order(abs(results$logFC), decreasing=T), ][1:5,"ID"]

#for (i in top_genes){
#	gene_name = subset(gene_description, ID == i)$Name
#	gene_levels <- data.frame("ID"=colnames(Vfinal$E), 
#		"BMI"=covariates[colnames(Vfinal$E), "log2BMI"],
#		"categories"=covariates[colnames(Vfinal$E), "BMI_categories"],
#		"TPM"=unlist(TPM[i,colnames(Vfinal$E)]),
#		"voom"=unlist(Vfinal$E[i,]),
#		"residuals"=as.numeric(unlist(full_model_residuals[i,colnames(Vfinal$E)]))
#		)
#	if(contrast_is_dicho){violin_plot(gene_name, gene_levels)
#	} else { scatter_plot(gene_name, gene_levels) }
#}
