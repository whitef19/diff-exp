#--------------------------------
# Differential expression analysis with SVA and limma-voom  
# by Frederique White

usage = function(errM) 
{
	cat("\nUsage : Rscript de-analysis.R -c <Value> -o <Value>\n")
	cat("       -c      : config\n")
	cat("       -o      : output directory\n")
	cat("       -h      : this help\n\n")
	stop(errM)
}

edger_norm_nomm <- function(counts) 
{
	
	dataObject = edgeR::DGEList(counts=counts)
	dataObject = edgeR::calcNormFactors(dataObject, method="TMM")
	
	libsize = dataObject$samples$lib.size
	normfactors = dataObject$samples$norm.factors
	libsize = libsize * normfactors
	
	tmm <- counts / libsize * 1e6
	cpm <- data.frame(edgeR::cpm(dataObject))
	colnames(cpm) <- colnames(counts)

	return(list(tmm=tmm, cpm=cpm, data_object=dataObject))
}

edger_norm <- function(counts, model_matrix) 
{
	
	dataObject = edgeR::DGEList(counts=counts)
	dataObject = edgeR::calcNormFactors(dataObject, method="TMM")
	
	libsize = dataObject$samples$lib.size
	normfactors = dataObject$samples$norm.factors
	libsize = libsize * normfactors
	
	tmm <- counts / libsize * 1e6
	cpm <- data.frame(edgeR::cpm(dataObject))
	colnames(cpm) <- colnames(counts)

	Voom <- limma::voom(dataObject, design=model_matrix)
	
	return(list(tmm=tmm, cpm=cpm, data_object=dataObject, voom=Voom))
}

how_much_svs = function(voom, design, full_model, full_model_matrix, null_model_matrix) {
	

	### estimate N SVs
	residuals <- t(resid(lm(as.formula(paste0('t(voom)', full_model)), data=design)))
	isvaResult = isva::EstDimRMT(residuals, F)
	numSVs <- isvaResult$dim + 1

	### try if the model is converging with numSVs SVs, if not, rerun with 1 SV less 
	res <- try(SmartSVA::smartsva.cpp(voom, full_model_matrix, mod0=null_model_matrix, n.sv=numSVs, alpha=1, B=200, VERBOSE=F))
	if (inherits(res, "try-error")) {
		while (inherits(res, "try-error")){
			numSVs = numSVs - 1
			if(numSVs < 2) { stop("SVA model is not converging") }
			res <- try(SmartSVA::smartsva.cpp(voom, full_model_matrix, mod0=null_model_matrix, n.sv=numSVs, alpha=1, B=200, VERBOSE=F))
		}
	}
	
	cat(paste0(" model is converging\n => estimation of ", numSVs," SVs\n"), file=log)

	return(numSVs)
}


correlation_matrix <- function(design, numSVs) {
	
<<<<<<< HEAD
	variable_of_interest = c("Patient_Age","gestage_numeric","Sex","GravidStatusDich","smoking",
		"BMI_V1","GDM_IADSPG","Matsuda","Stumvoll","PE","GH",
		"rRNA_rate","high_quality_reads","intergenic_reads","intragenic_reads","low_quality_reads","mapped_reads","median_avg_transcript_coverage")
=======
	variable_of_interest = c("Patient_Age","gestage_num","sex_num","primigravid","Maternal_smoking",
		"Visit_BMI_V1","GDM_IADPSG","Matsuda","Stumvoll","PE","HTAg",
		"Batch","RIN_Novogene","Mean_Quality")
>>>>>>> 02e35f0945fcdba204505603ab13d0e9c7b31f3d
	df = design[, c(variable_of_interest, paste0("SV",1:numSVs))]
	#df$sex_num <- as.numeric(df$sex_num)
	#df$primigravid <- as.numeric(df$primigravid)
	#df$Batch <- as.numeric(df$Batch)
	master <- cor(na.omit(df))
	pdf(paste0(out_path, "/sva_corrplot.pdf"), width=12, height=12)
	corrplot(master, tl.col="black",method = 'color')
<<<<<<< HEAD
	corrplot(master[1:18,19:ncol(master)], tl.col="black",method = 'color')
=======
	corrplot(master[1:length(variable_of_interest),(length(variable_of_interest)+1):ncol(master)], tl.col="black",method = 'color')
>>>>>>> 02e35f0945fcdba204505603ab13d0e9c7b31f3d
	dev.off()
}

do_a_pca <- function(df, design, config) {

	PCA <- prcomp(t(df))
	PCs <- PCA$x
	#PCs <- data.frame( cbind(ID=rownames(PCs), PCs, conditions=design[rownames(PCs),config$parameters$contrast], data.frame(t(df[, rownames(PCs)]))) )
	PCs <- data.frame( cbind(ID=rownames(PCs), PCs, data.frame(t(df[, rownames(PCs)]))) )
	#if (nrow(data.frame(table(PCs$conditions))) >= 5 ) { PCs$contrast = cut(PCs$conditions, quantile(PCs$conditions), include=T) } else { PCs$contrast = PCs$conditions}
	var <- summary(PCA)$importance[2,]
	Var <- data.frame("x"=1:length(var), "var"=as.vector(var)*100)
	return(list(PC=PCs, var=Var, loadings=data.frame(PCA$rotation)))
}

complete_pca <- function(PCA, pdf_file, normalization) {

	top_loadings <- rownames(PCA$loadings[order(abs(PCA$loadings$PC1), decreasing=T), ])[1:3]
	#melted_top_loadings <- reshape2::melt(PCA$PC[, c("ID","contrast", "PC1", "PC2", top_loadings)], id=c("ID","contrast", "PC1", "PC2"))
	melted_top_loadings <- reshape2::melt(PCA$PC[, c("ID","PC1", "PC2", top_loadings)], id=c("ID","PC1", "PC2"))

	pdf(pdf_file, width=8, height=6)
	scree <- ggplot(PCA$var[1:20,], aes(x=x, y=var)) + geom_bar(stat="identity", fill="black") + ggtitle(paste0("Scree plot (",normalization,")")) + xlab("PCs") + ylab("Variance (%)")
	#pca <- ggplot(PCA$PC, aes(x=PC1, y=PC2, color=contrast, label=rownames(PCA$PC))) + geom_point(alpha=0.1) + ggtitle(paste0("PCA (",normalization,")")) + xlab(paste0("PC1 (", PCA$var$var[1] ," %)")) + ylab(paste0("PC2 (", PCA$var$var[2] ," %)")) + theme(legend.position="none")
	pca <- ggplot(PCA$PC, aes(x=PC1, y=PC2, label=rownames(PCA$PC))) + geom_point(alpha=0.1) + ggtitle(paste0("PCA (",normalization,")")) + xlab(paste0("PC1 (", PCA$var$var[1] ," %)")) + ylab(paste0("PC2 (", PCA$var$var[2] ," %)")) + theme(legend.position="none")
	loads <- ggplot(PCA$loadings, aes(x=PC1, y=PC2, label=rownames(PCA$loadings))) + geom_text_repel() + geom_point(show.legend=F, alpha=0.2) + ggtitle(paste0("loadings (",normalization,")"))
	empty <- ggplot() + theme_void()
	a <- ggarrange(empty, scree, empty, ncol=3, widths=c(1,3,1))
	b <- ggarrange(pca, loads, ncol=2)
	print(ggarrange(a, b, ncol=1, heights=c(1,2)))
	#print(ggplot(melted_top_loadings, aes(x=contrast, y=value, group=contrast, fill=contrast)) + geom_boxplot() + geom_jitter() + facet_wrap(~variable, scale="free") + ggtitle("Top loadings") + theme(legend.position="bottom"))
	#print(ggplot(melted_top_loadings, aes(x=PC1, y=value, color=contrast)) + geom_point() + ggtitle("Top loadings") + facet_wrap(~variable, scale="free")+ theme(legend.position="bottom"))
	print(ggplot(melted_top_loadings, aes(x=PC1, y=value)) + geom_point() + ggtitle("Top loadings") + facet_wrap(~variable, scale="free")+ theme(legend.position="bottom"))
	dev.off()

}
