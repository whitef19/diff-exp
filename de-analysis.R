#--------------------------------
# Differential expression analysis with SVA and limma-voom  
# by Frederique White
# 2022-10-28 
# 
# Usage : Rscript de-analysis.R -c path_to_config_file -o output_directory
#
# 
# Libraries
suppressMessages(library(edgeR))
suppressMessages(library(limma))

#--------------------------------

usage = function(errM) {
	cat("\nUsage : Rscript de-analysis.R -c <Value> -o <Value>\n")
	cat("       -c      : config\n")
	cat("       -o      : output directory\n")
	cat("       -h      : this help\n\n")
	stop(errM)
}

open_config = function(config_path) {
	if (!(file.exists(config_path))) {
		usage("Error : config file not found")
	}
	
	config = read.csv(config_path, sep="\t", header=T, row.names=1)
	if(ncol(config) != 2) {
		stop("config file must have 3 columns: variable names, types and values")
	} 

	colnames(config) <- c("variable_type", "variable_value")
	paths = data.frame(t(subset(config, variable_type == "file", "variable_value")))
	parameters = data.frame(t(subset(config, variable_type == "parameter", "variable_value")))

	# verify paths
	for (file in paths) {
		if (!is.na(file)){
			if (!(file.exists(file))) {
				stop(paste0("Error : ",file," file not found"))
			}
		}
	}

	return(list(paths=paths, parameters=parameters))	
}


filtering = function(config, design, included_samples, normalized_count_included){
	
	## read count file
	count_df = read.csv(config$paths$read_count_matrix, sep="\t", check.names=F, row.names=1)
	cat(paste0("=> initial number of genes: ", nrow(count_df), "\n"), file=log)
	
	## remove first column that contains description
	count_df$Description <- NULL
	
	## Apply count thresholds
	cat(paste0("Gene filter thresholds:\n    counts: ", config$parameters$count_threshold,
						"\n    sample fraction:   ", config$parameters$sample_frac_threshold,"\n"), file=log)
	## to be included, genes must have at least count_threshold in at least sample_frac_threshold of samples
	passed_count_threshold <- apply(count_df, 1, function(x) (sum(x >= as.numeric(config$parameters$count_threshold))/length(x)) >= as.numeric(config$parameters$sample_frac_threshold))

	## Apply normalized count thresholds

	if ( normalized_count_included ) {
		cat(paste0("    normalized counts: ", config$parameters$normalized_count_threshold,"\n"), file=log)
		normalized_count_df = read.csv(config$paths$normalized_count_matrix, sep="\t", check.names=F, row.names=1)[, -1]
		if (length(intersect(colnames(count_df), colnames(normalized_count_df))) != ncol(count_df)) { stop("samples are different between the counts file and the tpm file") }
		passed_normalized_threshold <- apply(normalized_count_df, 1, function(x) (sum(x >= as.numeric(config$parameters$normalized_count_threshold))/length(x)) >= as.numeric(config$parameters$sample_frac_threshold))
		passed_genes <- passed_count_threshold & passed_normalized_threshold
	} else {
		passed_genes <- passed_count_threshold
	}

	count_df <- as.data.frame(count_df[passed_genes,])
	cat(paste0("=> passed genes: ", nrow(count_df),"\n"), file=log)

	included_genes <- rownames(count_df)

	# Apply mappability threshold
	if (file.exists(config$paths$mappability_scores)) {
		mappability_scores <- read.table(config$paths$mappability_scores, header=F, col.names=c("gene_id","score"))
		high_mappability_genes <- subset(mappability_scores, score >= config$parameters$mappability_threshold)$gene_id
		included_genes <- intersect(high_mappability_genes, included_genes)
		count_df <- count_df[included_genes, ]
		cat(paste0("# mappability score threshold: ", config$parameters$mappability_threshold,"\n=> high mappability genes: ", nrow(count_df),"\n"), file=log)
	}

	## sample filtering
	cat(paste0("=> initial number of samples: ", ncol(count_df), "\n"),file=log)

	## include only samples without any missing covariates
	design = design[included_samples, ]
	model <- paste0("~", config$parameters$model)
	mm <- model.matrix(eval(parse(text=model)), data=design) # model.matrix remove rows with missing values 

	included_samples <- rownames(mm)
	cat(paste0("=> number of included samples: ", length(included_samples),"\n"),file=log)
	
	count_df = count_df[, included_samples] 
	clean_count_df <- cbind("ID"=rownames(count_df), count_df)
	clean_file_name = paste0(out_path,"/clean_counts.tsv")
	write.table(clean_count_df, file=paste0(out_path,"/clean_counts.tsv"), sep="\t", quote=F, row.names=F)
	
	if ( normalized_count_included ) {
		return(list(count_df=count_df, included_samples=included_samples, included_genes=included_genes, normalized=normalized_count_df[included_genes, included_samples]))
	} else {
		return(list(count_df=count_df, included_samples=included_samples, included_genes=included_genes))
	}
}


run_sva = function(dataObject, design, config){

	suppressMessages(library(isva))
	suppressMessages(library(SmartSVA))

	model = paste0("~ ", config$parameters$model )
	mm <- model.matrix(eval(parse(text=model)), data=design) 
	cat(paste0("# Compute surrugate variables\n", " Model: ", model ,"\n"), file=log)
	
	V <- limma::voom(dataObject, design=mm)
	normalized_count_df <- V$E
	
	## estimate N SVs
	Yr <- t(resid(lm(as.formula(paste0('t(normalized_count_df)', model)), data=design)))
	isvaResult = isva::EstDimRMT(Yr, F)
	numSVs <- isvaResult$dim + 1
	cat(paste0(" trying with ", numSVs," SVs\n"), file=log)

	## try if the model is converging with numSVs SVs, if not, rerun with 1 SV less 
	res <- try(SmartSVA::smartsva.cpp(normalized_count_df, mm, mod0=NULL, n.sv=numSVs, alpha=1, B=200, VERBOSE=F))
	if (inherits(res, "try-error")) {
		while (inherits(res, "try-error")){
			numSVs = numSVs - 1
			if(numSVs < 30) { stop("SVA model is not converging") }
			res <- try(SmartSVA::smartsva.cpp(normalized_count_df, mm, mod0=NULL, n.sv=numSVs, alpha=1, B=200, VERBOSE=F))
		}
	}
	
	## run sva
	cat(paste0(" model is converging\n => estimation of ", numSVs," SVs\n"), file=log)
	SVObject <- SmartSVA::smartsva.cpp(normalized_count_df, mm, mod0=NULL, n.sv=numSVs, alpha=1, B=200, VERBOSE=F)

	## concatenate SVs to covariate dataframe 
	SVs <- as.data.frame(SVObject$sv)
	names(SVs) <- paste0("SV",1:ncol(SVs))
	design <- cbind(design, SVs)
	design_to_write <- cbind(rownames(design), design)
	write.table(design, file=paste0(out_path,"/design_with_SVs.tsv"))
	
	full_model = paste0("~", model, "+", paste0("SV",1:numSVs, collapse="+"))
	return(list(design=design, model=full_model))
}


run_limma = function(dataObject, design, config, normalized_count_included) {

	## create covariate matrix to input in limma
	mm = model.matrix(eval(parse(text=full_model)), data=design)
	
	## rerun voom with full model
	Vfinal = limma::voom(dataObject, design=mm)
	voom_count_df = Vfinal$E
	write.table(cbind("ID"=rownames(voom_count_df), voom_count_df), file=paste0(out_path, "/voom_normalized_counts.tsv" ), sep="\t", quote=F, row.names=F)

	## run the regressions
	fit <- limma::lmFit(Vfinal, mm)
	## run eBayes() to borrow information across genes
	ebayesFit <- limma::eBayes(fit)
	## Extract limma results
	results <- limma::topTable(ebayesFit, coef=config$parameters$contrast, number=Inf, sort.by="p")
	results <- cbind(ID=rownames(results), results)

	## Annotate limma results
	results$ID <- rownames(results)
	if (! is.na(config$paths$gene_annotations)) {
		gene_description = read.csv(config$paths$gene_annotations, sep="\t")
		results <- merge(gene_description, results, by="ID", all.y=T)
	}

	if (normalized_count_included){
		results <- merge(results, data.frame("ID"=rownames(normalized_df), "avgTPM"=rowMeans(normalized_df)), by="ID", all.x=T)
	}

	results <- results[order(results$P.Value),]
	write.table(results, file=paste0(out_path,"/limma_results.tsv"), sep="\t", quote=F, row.names=F)
	return(results)
}


create_figures = function(results, normalized_count_included){
	suppressMessages(library(ggplot2))
	suppressMessages(require(ggrepel))
	
	results$absFC = abs(results$logFC)
	logFC_threshold <- (mean(results$absFC) + 5 * sd(results$absFC)) 
	top_hits <- subset(results, absFC > logFC_threshold)
	results$color <- ifelse(results$absFC > logFC_threshold, "#22908C", "#F2BB05")
	results$color <- ifelse(results$P.Value > 0.05, "black", results$color)

	cat(paste0("logFC threshold mean + 5 SDs: ",logFC_threshold,  "\n",
					"==> number of top hits: ",nrow(top_hits),"\n"), file=log)
	write.table(top_hits, file=paste0(out_path,"/top_limma_results.tsv"), sep="\t", quote=F, row.names=F)

	mean_logfc <- mean(results$absFC)
	sd_logfc <- sd(results$absFC)
	
	pdf(paste0(out_path, "/figures.pdf"), width=12, height=6)

	#print(ggplot(results, aes(x=ID, y=absFC)) + geom_point(color="#22908C", alpha=0.5) + ylab("absolute log2 fold change") +
	#	geom_hline(yintercept=mean_logfc, color="#22908C",linetype="dashed") +
	#	geom_hline(yintercept=c(sapply(1:5,function(x) c( mean_logfc + x * sd_logfc )), color="#F2BB05" ,linetype="dashed") + 
	#	annotate(geom="text", x=nrow(results)/2, y=mean_logfc + 6*sd_logfc, label=paste0("avg. log2 fold change + 5*SD: ",signif(logFC_threshold, digits=3)), color="black")))
	# volcano plot	
	volcano <- ggplot(results, aes(x=logFC, y=-log10(P.Value), color=color)) + geom_point(show.legend=F) + 
		ylab("-log10 nominal p-value") + xlab("log2 fold change") + 
		xlim(-max(results$absFC), max(results$absFC)) +
		ggtitle(paste0("~ ",config$parameters$model,"\n + ",config$parameters$number_surrogate_variables, " SVs","\n n = ",length(included_samples)," samples")) + 
		geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), color="grey", linetype="dashed") +
		geom_hline(yintercept=c(-log10(0.05), -log10(0.001)), color="grey",linetype="dashed") +
		theme(plot.title=element_text(size=7)) + scale_color_identity()

	if (nrow(subset(results, abs(logFC) > logFC_threshold & P.Value < 0.001)) > 0){
		if ( ! is.na(config$paths$gene_annotations)) {
			print(volcano + geom_text_repel(data =subset(results, abs(logFC) > logFC_threshold & P.Value < 0.001), aes(label=Name, size=1.5), show.legend=F))
		} else {
			print(volcano + geom_text_repel(data=subset(results, abs(logFC) > logFC_threshold & P.Value < 0.001), aes(label=ID, size=1.5), show.legend=F))
		}
	} else {print(volcano)}

	dev.off()

}





### main 

ARG = commandArgs(trailingOnly = T)
if(length(ARG)<2) { usage("missing arguments")}

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
config = open_config(config_path)

## read design file ans list of sample to included in analysis
design = read.csv(config$paths$design, sep="\t", check.names=F, row.names=1)
included_samples = read.csv(config$paths$included_samples, header=F)$V1

## determine if normalized counts are included
normalized_count_included = ifelse(is.na(config$paths$normalized_count_matrix), FALSE, TRUE)

## which count file to use
if (! is.null(config$paths$clean_counts)){
	
	clean_count_file = 	config$paths$clean_counts

} else if (file.exists(paste0(out_path, "/clean_counts.tsv" ))) {

	cat(paste0(out_path, "/clean_counts.tsv"," exist ! Do you want to use it (y) or replace it (n): "));
	answer = readLines("stdin",n=1);

	if (answer == "y") { clean_count_file = paste0(out_path, "/clean_counts.tsv")}
	else if (answer == "n") { clean_count_file = NA } 
	else {stop()}
} else {
	clean_count_file = NA
}


## create or read clean count file
if (is.na(clean_count_file)) {
	clean = filtering(config, design, included_samples, normalized_count_included)
	count_df = clean$count_df
	included_samples = clean$included_samples
	included_genes = clean$included_genes
	if (normalized_count_included) {normalized_df = clean$normalized }

} else {
	cat(paste0("using ", clean_count_file, "\n"), file=log)
	count_df = read.csv(clean_count_file, sep="\t", row.names=1, check.names=F)
	included_genes = rownames(count_df)
	included_samples = intersect(included_samples, colnames(count_df))
	count_df = count_df[, included_samples]
	if(normalized_count_included){
		normalized_df = read.csv(config$paths$normalized_count_matrix, sep="\t", check.names=F, row.names=1)
		normalized_df = normalized_df[included_genes, included_samples]
	}
}


## filter design
design = design[included_samples, ]

## normalize counts
dataObject = edgeR::DGEList(counts=count_df)
dataObject = edgeR::calcNormFactors(dataObject)


## run surrogate variable analysis
if (as.logical(config$parameters$run_SVA)) {
	from_sva = run_sva(dataObject, design, config)
	design = from_sva$design
	full_model = from_sva$model
} else if (config$parameters$number_surrogate_variables > 0) {
	full_model = paste0("~", config$parameters$model, "+", paste0("SV",1:config$parameters$number_surrogate_variables, collapse="+"))
} else {
	full_model = paste0("~", config$parameters$model)
}
cat(paste0("# Full Model: ", full_model ,"\n"), file=log)

## Differential expression analysis 
results = run_limma(dataObject, design, config, normalized_count_included)

create_figures(results, normalized_count_included)




violin_plot <- function(gene_name, gene_levels){
	suppressMessages(library(ggpubr))
	p1 <- ggplot(gene_levels, aes(x=categories, y=normalized, fill=categories)) + xlab("") + ylab(paste0("TPM ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	p2 <- ggplot(gene_levels, aes(x=categories, y=voom, fill=categories)) + xlab("") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	p3 <- ggplot(gene_levels, aes(x=categories, y=residuals, fill=categories)) + xlab("") + ylab(paste0("Full model residuals ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	print(ggarrange(p1, p2, p3, ncol=3))
}

scatter_plot <- function(gene_name, gene_levels){
	suppressMessages(library(ggpubr))
	p1 <- ggplot(gene_levels, aes(x=BMI, y=TPM)) + xlab("log2 BMI") + ylab(paste0("TPM ", gene_name )) + geom_point(alpha=0.5,color="#22908C") 
	p2 <- ggplot(gene_levels, aes(x=BMI, y=voom)) + xlab("log2 BMI") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
	p3 <- ggplot(gene_levels, aes(x=BMI, y=residuals)) + xlab("log2 BMI") + ylab(paste0("Full model residuals ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
	print(ggarrange(p1, p2, p3, ncol=3))
}

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
