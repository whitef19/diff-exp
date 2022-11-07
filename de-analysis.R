#--------------------------------
# Differential expression analysis with SVA and limma-voom  
# by Frederique White
# 2022-10-28 
# 
# Usage : Rscript de-analysis.R -c path_to_config_file -o output_directory
#
# 
# Libraries
#suppressMessages(library(edgeR))
#suppressMessages(library(isva))
#suppressMessages(library(SmartSVA))
#suppressMessages(library(limma))
#suppressMessages(library(ggplot2))
#suppressMessages(library(ggpubr))
#suppressMessages(require(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
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
		if (!(file.exists(file))) {
			stop(paste0("Error : ",file," file not found"))
		}
	}

	return(list(paths=paths, parameters=parameters))	
}


filtering = function(config, design, included_samples){
	
	count_df = read.csv(config$paths$read_count_matrix, sep="\t", check.names=F, row.names=1)[,-1]
	tpm_df = read.csv(config$paths$tpm_count_matrix, sep="\t", check.names=F, row.names=1)[, -1]
	
	if (length(intersect(colnames(count_df), colnames(tpm_df))) != ncol(count_df)) { 
		stop("samples are different between the counts file and the tpm file") 
	}
	

	## gene filtering
	cat(paste0("=> initial number of genes: ", nrow(count_df), "\n"), file=log)
	cat(paste0("# Gene filter thresholds:\n    counts: ", config$parameters$count_threshold,
									"\n    TPM: ", config$parameters$tpm_threshold,
									"\n    sample fraction:   ", config$parameters$sample_frac_threshold,"\n"), file=log)
	## Apply counts and transcript per million (TPM) thresholds
	# x takes place of each rows
	passed_count_threshold <- apply(count_df, 1, function(x) sum(x >= config$parameters$count_threshold)/length(x) >= config$parameters$sample_frac_threshold)
	passed_tpm_threshold <- apply(tpm_df, 1, function(x) sum(x >= config$parameters$tpm_threshold)/length(x) >= config$parameters$sample_frac_threshold)
	passed_genes <- passed_count_threshold & passed_tpm_threshold
	count_df <- as.data.frame(count_df[passed_genes,])
	cat(paste0("=> passed genes: ", nrow(count_df),"\n"), file=log)

	included_genes <- rownames(count_df)

	# Apply mappability threshold
	if (file.exists(config$paths$mappability_scores)) {
		mappability_scores <- read.table(config$paths$mappability_scores, header=F, col.names=c("gene_id","score"))
		high_mappability_genes <- subset(mappability_scores, score >= config$parameters$mappability_threshold)$gene_id
		included_genes <- intersect(high_mappability_genes, included_genes)
		count_df <- subset(count_df, rownames(count_df) %in% high_mappability_genes)
		cat(paste0("# mappability score threshold: ", config$parameters$mappability_threshold,"=> high mappability genes: ", nrow(count_df),"\n"), file=log)
	}

	## sample filtering
	cat(paste0("=> initial number of samples: ", ncol(count_df), "\n"),file=log)

	## include only samples without any missing covariates
	design = design[included_samples, ]
	model <- paste0("~", config$parameters$model)
	mm <- model.matrix(eval(parse(text=model)), data=design) # model.matrix remove rows with missing values 

	print("ok")

	included_samples <- rownames(mm)
	cat(paste0("=> number of included samples: ", length(included_samples),"\n"),file=log)
	
	count_df = count_df[, included_samples] 
	clean_count_df <- cbind("ID"=rownames(count_df), count_df)
	clean_file_name = paste0(out_path,"/clean_counts.tsv")
	print("ok...")

	if ((file.exists(clean_file_name))) {usage("Error : config file not found")}
	write.table(clean_count_df, file=paste0(out_path,"/clean_counts.tsv"), sep="\t", quote=F, row.names=F)

	return(list(count_df=count_df, included_samples=included_samples, included_genes=included_genes))
}


SVA = function(dataObject, design, config){

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

if(! dir.exists(out_path)) { dir.create(out_path) }
log = file(paste0(out_path,"/errors.o"), open="a")

config <- open_config(config_path)

design = read.csv(config$paths$design, sep="\t", check.names=F, row.names=1)
included_samples = read.csv(config$paths$included_samples, header=F)$V1
gene_description = fread(config$paths$read_count_matrix, sep="\t", select=c("Name", "Description"))
colnames(gene_description) = c("ID", "Name")

## gene and sample filtering
## look if clean file already exist 
answer = "n"
clean_count_file = paste0( out_path, "/clean_counts.tsv" )
if (file.exists(clean_count_file)) {
	cat(paste0(clean_count_file," exist ! Do you want to use it (y) or replace it (n): "));
	answer = readLines("stdin",n=1);
}

if( answer == "n") {
	clean = filtering(config, design, included_samples)
	count_df = clean$count_df
	included_samples = clean$included_samples
	included_genes = clean$included_genes
} else if (answer == "y"){
	count_df = read.csv(clean_count_file, sep="\t", row.names=1)
	included_genes = rownames(count_df)
	included_samples = intersect(included_samples, colnames(count_df))
	count_df = count_df[, included_samples]
} else { stop() }

print("that is done")

## filter design
design = design[included_samples, ]

## normalize counts
dataObject = edgeR::DGEList(counts=counts)
dataObject = edgeR::calcNormFactors(dataObject)


## run surrogate variable analysis
if (as.boolean(config$parameters$run_SVA)) {
	from_sva = SVA(dataObject, design, config)
	design = from_sva$design
	full_model = from_sva$model
} else if (config$parameters$numSVs > 0) {
	full_model = paste0("~",config$parameters$model, "+", paste0("SV",1:numSVs, collapse="+"))
} else {
	full_model = paste0("~ ", config$parameters$model )
}
cat(paste0("# Full Model: ", full_model ,"\n"), file=log)


stop("you've reach the end")


## Differential expression analysis 
mm <- model.matrix(eval(parse(text=full_model)), data=design)
Vfinal <- limma::voom(dataObject, design=mm)
normalized_count_df <- Vfinal$E
write.table(cbind("ID"=rownames(normalized_count_df), normalized_count_df), file=paste0(out_path, "/voom_normalized_counts.tsv" ), sep="\t", quote=F, row.names=F)

## run the regressions
fit <- limma::lmFit(Vfinal, mm)
## run eBayes() to borrow information across genes
ebayesFit <- limma::eBayes(fit)
## Extract limma results
results <- limma::topTable(ebayesFit, coef=config$parameters$contrast, number=Inf, sort.by="p")

## Annotate limma results
results$ID <- rownames(results)
results <- merge(gene_description, results, by="ID", all.y=T)
results <- merge(results, data.frame("ID"=rownames(TPM), "avgTPM"=rowMeans(TPM)), by="ID", all.x=T)
results <- results[order(results$P.Value),]

write.table(results, file=results_out, sep="\t", quote=F, row.names=F)
if (included_sex !="both") {
	write.table(results, file=paste0(contrast,".",dataset,".",samples_subset,included_sex,".full_table.tsv" ), sep="\t", quote=F, row.names=F)
} else {
	write.table(results, file=paste0(contrast,".",dataset,".",samples_subset,"full_table.tsv" ), sep="\t", quote=F, row.names=F)}




#----------------------------------
# Post-processing

determine_threshold <- function(x) { return((mean(x) + 5 * sd(x)))}

logFC_threshold <- determine_threshold(abs(results$logFC))
top_hits <- subset(results, abs(logFC) > logFC_threshold)

cat(paste0("# logFC threshold mean + 5 SDs: ",logFC_threshold,  "\n",
				"==> number of top hits: ",nrow(top_hits),"\n"),
				file=log
				)

write.table(top_hits, file=top_hits_out, sep="\t", quote=F, row.names=F)
if (included_sex !="both") {
	write.table(top_hits, file=paste0(contrast,".",dataset,".",samples_subset, included_sex,".top_hits.tsv" ), sep="\t", quote=F, row.names=F)
} else {
	write.table(top_hits, file=paste0(contrast,".",dataset,".",samples_subset,"top_hits.tsv" ), sep="\t", quote=F, row.names=F)}




# visualize results

pdf(plot_out, width=12, height=6)

# distribution
ggplot(results, aes(x=P.Value)) + geom_histogram(bins=50, fill="#22908C") + xlab("nominal p-value")
ggplot(results, aes(x=logFC)) + geom_histogram(bins=50,fill="#22908C") + xlab("log2 FC")

mean_logFC <- mean(abs(results$logFC))
sd_logFC <- sd(abs(results$logFC))
ggplot(results, aes(x=Name, y=abs(logFC))) + geom_point(color="#22908C", alpha=0.5) + ylab("absolute log2 FC") +
	geom_hline(yintercept=mean_logFC, color="#22908C",linetype="dashed") +
	geom_hline(yintercept=c(sapply(1:5,function(x) c( mean_logFC + x * sd_logFC )), color="#F2BB05" ,linetype="dashed") + 
	annotate(geom="text", x=nrow(results)/2, y=mean_logFC+6*sd_logFC, label=paste0("avg. log2FC+5*SD: ",signif(logFC_threshold, digits=3)), color="black"))

# volcano plot
results$color <- ifelse(abs(results$logFC) > logFC_threshold, "#22908C", "#F2BB05")
ggplot(results, aes(x=logFC, y=-log10(P.Value), color=color)) + geom_point(show.legend=F) + 
	ylab("-log10 nominal p-value") + xlab("log2 FC") + 
	xlim(-max(abs(results$logFC)), max(abs(results$logFC))) +
	ggtitle(paste0("~ ",contrast, " + ", null_model,"\n + ",numSVs, " SVs","\n n = ",length(sample_list)," samples\n",dataset, " ",samples_subset,"\n","sex included:",included_sex)) + 
	geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), color="grey", linetype="dashed") +
	geom_hline(yintercept=c(-log10(0.05), -log10(0.001)), color="grey",linetype="dashed") +
	geom_text_repel(data =subset(results, abs(logFC) > logFC_threshold & P.Value < 0.001),aes(label=Name, size=1.5), show.legend=F) + 
	theme(plot.title=element_text(size=7))+scale_color_identity()



violin_plot <- function(gene_name, gene_levels){
	p1 <- ggplot(gene_levels, aes(x=categories, y=TPM, fill=categories)) + xlab("") + ylab(paste0("TPM ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	p2 <- ggplot(gene_levels, aes(x=categories, y=voom, fill=categories)) + xlab("") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	p3 <- ggplot(gene_levels, aes(x=categories, y=residuals, fill=categories)) + xlab("") + ylab(paste0("Full model residuals ", gene_name )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed")
	print(ggarrange(p1, p2, p3, ncol=3))
}

scatter_plot <- function(gene_name, gene_levels){
	p1 <- ggplot(gene_levels, aes(x=BMI, y=TPM)) + xlab("log2 BMI") + ylab(paste0("TPM ", gene_name )) + geom_point(alpha=0.5,color="#22908C") 
	p2 <- ggplot(gene_levels, aes(x=BMI, y=voom)) + xlab("log2 BMI") + ylab(paste0("Voom-normalized expression ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
	p3 <- ggplot(gene_levels, aes(x=BMI, y=residuals)) + xlab("log2 BMI") + ylab(paste0("Full model residuals ", gene_name )) + geom_point(alpha=0.5, color="#22908C") 
	print(ggarrange(p1, p2, p3, ncol=3))
}

contrast_is_dicho <- nrow(as.data.frame(table(covariates[, contrast]))) == 2
top_genes <- results[order(abs(results$logFC), decreasing=T), ][1:5,"ID"]

for (i in top_genes){
	gene_name = subset(gene_description, ID == i)$Name
	gene_levels <- data.frame("ID"=colnames(Vfinal$E), 
		"BMI"=covariates[colnames(Vfinal$E), "log2BMI"],
		"categories"=covariates[colnames(Vfinal$E), "BMI_categories"],
		"TPM"=unlist(TPM[i,colnames(Vfinal$E)]),
		"voom"=unlist(Vfinal$E[i,]),
		"residuals"=as.numeric(unlist(full_model_residuals[i,colnames(Vfinal$E)]))
		)
	if(contrast_is_dicho){violin_plot(gene_name, gene_levels)
	} else { scatter_plot(gene_name, gene_levels) }
}

dev.off()


#fc_threshold <- mean(fc_distribution$slopes) + 5 * sd(fc_distribution$slopes)
#print(fc_threshold)
#x = fc_distribution$slopes
#is_outlier <- function(x) { return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))}
#results$colors <- ifelse(results$logFC > fc_threshold, "darkcyan", ifelse(results$logFC < -fc_threshold, "deeppink4", "grey") )
#results$Colors <- ifelse(results$P.Value < 0.05, results$colors, "grey")
#results$colors <- ifelse(results$logFC > fc_threshold, "black", ifelse(results$logFC < -fc_threshold, "black", "grey") )
#results$Colors <- ifelse(results$P.Value < 0.001, results$colors, "grey")
#ggplot(results, aes(x=abs(logFC), y=B)) + geom_point(size=2, alpha=0.3)
#contrast_df <- data.frame("ID"=rownames(covariates), "contrast"=covariates[, contrast])
#contrast_is_dicho <- nrow(as.data.frame(table(contrast_df$contrast))) == 2
if(2 > 3){
	print("Dichotomous analysis")

	# volcano plot
	print(ggplot(results, aes(x=logFC, y=-log10(P.Value), color=Colors)) + geom_point(show.legend=F) + xlab("log2 fold change") + ylab("-log10(p-value)") + xlim(-max(abs(results$logFC)), max(abs(results$logFC))) + 
	  ggtitle(paste0(model,"\n",dataset, " ",samples_subset,"\n","sex:",included_sex)) + scale_color_identity() + theme(plot.title=element_text(size=7)) + geom_vline(xintercept=c(-fc_threshold, fc_threshold), color="grey", linetype="dashed") +
	  geom_hline(yintercept=-log10(0.001), color="grey",linetype="dashed") +
	  geom_text_repel(data =subset(results, abs(logFC) > fc_threshold & P.Value < 0.001),aes(label=Name, size=2), show.legend=F))

	print(ggplot(fc_distribution, aes(x=genes, y=slopes)) + geom_point(size=2, alpha=0.15) + ylab("absolute log2 fold change") + geom_hline(yintercept=c(sapply(1:5,function(x) c(mean(fc_distribution$slopes)+x*sd(fc_distribution$slopes)) )), linetype="dashed") + 
		geom_hline(yintercept=mean(fc_distribution$slopes), color="red",linetype="dashed") +
		annotate(geom="text", x=2000, y=mean(fc_distribution$slopes)+6*sd(fc_distribution$slopes), label=paste0("avg. log2FC+5SD: ",signif(fc_threshold, digits=3)), color="black"))

	contrast_df$group <- ifelse(contrast_df$contrast == 0, "Normal", "Obese") 
	num <- as.data.frame(table(contrast_df$group))
	colnames(num) <- c("group","num")
	num$name <- paste0(num$group,"\nn=",num$num)
	contrast_df <- merge(contrast_df, num, by="group")

	library(reshape2)
	top_genes <- results[1:10, "ID"]
	expr <- data.frame(t(Vfinal$E[top_genes, ]))
	expr$ID <- rownames(expr)
	expr <- merge(expr, contrast_df, by="ID")
	melted <- melt(expr, id=c("ID", "contrast", "num", "name", "group"))	
	melted <- merge(melted, gene_description, by.x="variable", by.y="ID")
	print(ggplot(melted, aes(x=Name, y=value, fill=group)) + geom_boxplot() + xlab("Genes") +ylab("Voom-normalized expression")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

	residuals <- data.frame(t(full_model_residuals[top_genes, ]))
	residuals$ID <- rownames(residuals)

	residuals <- merge(residuals, contrast_df, by="ID")

	melted <- melt(residuals, id=c("ID", "contrast", "num", "name", "group"))	
	print("ok")

	melted <- merge(melted, gene_description, by.x="variable", by.y="ID")
	print(ggplot(melted, aes(x=Name, y=value, fill=group)) + geom_boxplot() + xlab("Genes") +ylab("Voom-normalized residuals from full model")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	print(Vfinal$E[1:5,1:5])
	print("ok")

	for (i in 1:5){
		top_gene <- results[i, "ID"]

		top_gene_expression <- data.frame("ID"=colnames(Vfinal$E), "expr"=unlist(Vfinal$E[top_gene,]), "residuals"=unlist(full_model_residuals[top_gene,]))
		tmp <- merge(top_gene_expression, contrast_df, by="ID")
		print(ggplot(tmp, aes(x=name, y=expr, fill=group)) + xlab("")+ylab(paste0("Voom-normalized expression ", results[i, "Name"] )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed"))
		print(ggplot(tmp, aes(x=name, y=residuals, fill=group)) + xlab("")+ylab(paste0("Voom-normalized residuals from full model ", results[i, "Name"] )) + geom_violin() + geom_boxplot(width=0.2, alpha=0.2)+stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .75, linetype = "dashed"))
	}

	dev.off()

}
if (3>4) {
	print("Continuous analysis")

	# volcano plot

	print(ggplot(results, aes(x=logFC, y=-log10(P.Value), color=Colors)) + geom_point(show.legend=F) + xlab("slope") + ylab("-log10(p-value)") + xlim(-max(abs(results$logFC)), max(abs(results$logFC))) +
	  ggtitle(paste0(model,"\n",dataset, " ",samples_subset,"\n","sex:",included_sex)) + scale_color_identity() + theme(plot.title=element_text(size=7)) + geom_vline(xintercept=c(-fc_threshold, fc_threshold), color="grey", linetype="dashed") +
	  geom_hline(yintercept=-log10(0.001), color="grey",linetype="dashed") +
	  geom_text_repel(data =subset(results, abs(logFC) > fc_threshold & P.Value < 0.001),aes(label=Name, size=2), show.legend=F))

	# threholds
	print(ggplot(fc_distribution, aes(x=genes, y=slopes)) + geom_point(size=2, alpha=0.15) + ylab("absolute slopes") + geom_hline(yintercept=c(sapply(1:5,function(x) c(mean(fc_distribution$slopes)+x*sd(fc_distribution$slopes)) )), linetype="dashed") + 
	  geom_hline(yintercept=mean(fc_distribution$slopes), color="red",linetype="dashed") +
	  annotate(geom="text", x=2000, y=mean(fc_distribution$slopes)+6*sd(fc_distribution$slopes), label=paste0("avg. slope+5SD: ",signif(fc_threshold, digits=3)), color="black"))
	
	print(ggplot(results, aes(x=order, y=scaled_slope)) + geom_point(size=2, alpha=0.15) + ylab("Z-score slopes") + geom_text_repel(data=subset(results, abs(scaled_slope) > 10), aes(label=Name)))
	
	for (i in 1:5){
		top_gene <- results[i, "ID"]
		top_gene_expression <- data.frame("ID"=colnames(Vfinal$E), "expr"=Vfinal$E[top_gene,])
		top_gene_residuals <- data.frame("ID"=colnames(full_model_residuals), "residuals"=full_model_residuals[top_gene,])

		top_gene_expression <- merge(top_gene_expression, contrast_df, by="ID")
		top_gene_expression <- merge(top_gene_expression, top_gene_residuals, by="ID")
		slope <- coef(fit)[top_gene, contrast]
		print(ggplot(top_gene_expression, aes(x=contrast, y=expr)) + xlab(contrast) + ylab(paste0("Voom-normalized expression ", results[i, "Name"] )) + geom_point(alpha=0.2) + annotate(geom="text", x=quantile(top_gene_expression$contrast,0.05), y=quantile(top_gene_expression$expr, 0.95), label=paste0("slope: ",signif(slope, digits=3)), color="black")) 
		print(ggplot(top_gene_expression, aes(x=contrast, y=residuals)) + xlab(contrast) + ylab(paste0("Voom residuals from full model ", results[i, "Name"] )) + geom_point(alpha=0.2) + annotate(geom="text", x=quantile(top_gene_expression$contrast,0.05), y=quantile(top_gene_expression$residuals, 0.95), label=paste0("slope: ",signif(slope, digits=3)), color="black")) 
	}
	dev.off()
}	


