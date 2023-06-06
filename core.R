#--------------------------------
# # by Frederique White

usage = function(errM) {
	cat("\nUsage : Rscript {de|sv}-analysis.R -d <Value> -c <Value> -o <Value>\n")
	cat("       -d      : design\n")
	cat("       -c      : config\n")
	cat("       -o      : output directory\n")
	cat("       -h      : this help\n\n")
	stop(errM)
}

read_command_line_arguments = function(){
	## default arg values
	output_path <<- "./"
	## get arg variables
	ARG = commandArgs(trailingOnly = T)
	if (length(ARG) < 4) { usage("missing arguments") }
	
	for (i in 1:length(ARG)) {
		if (ARG[i] == "-c") {
			config_path <<- ARG[i+1]
		} else if (ARG[i] == "-d") {
			design_path <<- ARG[i+1]
		} else if (ARG[i] == "-o") {
			output_path <<- ARG[i+1]
		} else if (ARG[i] == "-h") {
			usage("")
		}}
	if (! dir.exists(output_path)) { dir.create(output_path) }
	
	log_file <<- file(paste0(output_path,"/log.o"), open="a")
	cat(paste0("### Run started time: ",Sys.Date(), " at ", Sys.time(),"\n"), file=log_file)
	cat(paste0("### Design file used: ",design_path,"\n"), file=log_file)
}

read_config_file = function() {

	if (!(file.exists(config_path))) {
		usage("Error : config file not found")
	}
	
	get_extension <- unlist(strsplit(config_path, ".", fixed=T))
	extension <- get_extension[length(get_extension)]

	# config object is read only, set as global variable
	if ( extension == "csv" ) {
		config_df = read.csv(config_path, header=F, sep="\t")
		write.table(config_df, file=paste0(output_path, "/config_file_used.csv", sep="\t",quote=F,row.names=F))
		config <<- list()
		for (i in 1:nrow(config_df)){
			config[[i]] = unlist(config_df[i, "V2"])
		}
		names(config) = config_df$V1
		print(config)

	} else if (extension == "json") {
		config <<- fromJSON(file=config_path)	
		write(toJSON(config), paste0(output_path, "/config_file_used.json"))
	}

	# verify paths
	file_list = c("included_samples", "read_count_matrix", "normalized_count_matrix", "gene_annotations","mappability_scores")
	for (file in file_list) {
		file = unlist(config[file])
		if (!is.na(file)) {
			if (!(file.exists(file))) {
				stop(paste0("Error : ",file," file not found"))}}}

	cat(paste0("### see ", output_path, "/config_file_used.",extension," for files and parameters used for this run.\n"), file=log_file)
}

read_quantification_matrix = function(filename) {
	count_df = read.csv(filename, sep="\t", check.names=F, row.names=1)
	#cat(paste0("=> initial number of genes in ",basename(filename)," : ", nrow(count_df), "\n"), file=log_file)	
	## remove second column that contains description
	if (!is.numeric(count_df[1,1])){
		count_df[, 1] <- NULL
	}
	
	return(count_df)
}

gene_filtering = function(count_df, normalized_count_df){
	### Apply count thresholds
	passed_genes <- apply(count_df, 1, function(x) (sum(x >= as.numeric(config$count_threshold))/length(x)) >= as.numeric(config$sample_frac_threshold))
	if (exists("normalized_count_df", mode="list")) {
		passed_norm_thres <- apply(normalized_count_df, 1, function(x) (sum(x >= as.numeric(config$normalized_count_threshold))/length(x)) >= as.numeric(config$sample_frac_threshold))
		passed_genes <- (passed_genes & passed_norm_thres)
	}
	included_gene_list <- rownames(count_df)[passed_genes]
	cat(paste0("=> gene passing count thresholds: ", length(included_gene_list),"\n"), file=log_file)

	### Apply mappability threshold
	if (!is.na(config$mappability_scores)) {		
		mappability_scores <- read.table(config$mappability_scores, header=F, col.names=c("gene_id","score"))
		high_mappability_genes <- subset(mappability_scores, score >= config$mappability_threshold)$gene_id
		included_gene_list <- intersect(high_mappability_genes, included_gene_list)
		cat(paste0("=> gene passing mappability scores : ", length(included_gene_list),"\n"), file=log_file)
	}
	return(included_gene_list)
}

prepare_dataframes <- function(mode="dea"){

	## read quantification matrix
	print(config$read_count_matrix)
	count_df = read_quantification_matrix(config$read_count_matrix)
	
	## open normalized count file if provided
	if (!is.na(config$normalized_count_matrix)) {
		normalized_count_df = read_quantification_matrix(config$normalized_count_matrix)
	} else {
		normalized_count_df = NA
	} 
	
	## filter quantification matrix
	included_gene_list = gene_filtering(count_df, normalized_count_df)

	## get initial list of samples to include in analysis
	included_sample_df <<- read.csv(config$included_samples, sep="\t", header=T)

	## overlap with quantification matrix
	## update number of sample included
	included_sample_df <<- subset(included_sample_df, as.character(quantification_table_id) %in% colnames(count_df))
	

	## get list of covariates to include in analysis 
	full_model_list <<- unlist(strsplit(gsub(" ", "",config$full_model),"+", fixed=T))
	
	## add surrogate variable if any included
	if ((as.numeric(config$number_estimated_variables_to_include) > 0) & mode=="dea"){
		full_model_list <<- c(full_model_list, paste0(config$estimated_variables_method,1:as.numeric(config$number_estimated_variables_to_include)))
	}

	## read covariates dataframe and keep included samples and included variables
	design_df = read.csv(design_path, sep="\t", check.names=F, row.names=1)
	design_df = design_df[as.character(included_sample_df$design_id), full_model_list]
	## look for missing data in design dataframe
	if (any(is.na(design_df))) {
		missing_values = setdiff(rownames(design_df), rownames(na.omit(design_df)))
		cat(paste0("!!! samples with missing values in design: ", paste0(missing_values, collapse=","),"\n"), file=log_file)
		design_id_to_include = intersect(rownames(design_df), rownames(na.omit(design_df)))
	} else {
		design_id_to_include = rownames(design_df)
	}
	
	## overlap with design
	## update number of sample included
	included_sample_df <<- subset(included_sample_df, as.character(design_id) %in% design_id_to_include)
	cat(paste0("=> number of included samples: ", nrow(included_sample_df),"\n"),file=log_file)

	## update design and quantification matrix with final list of included samples and genes
	clean_count_df = count_df[included_gene_list, as.character(included_sample_df$quantification_table_id)]
	clean_design_df = design_df[as.character(included_sample_df$design_id), ]
	
	return (list(design_df=clean_design_df, count_df=clean_count_df))
}

edgeR_normalization <- function(count_df, design_df) {
	
	model_matrix = model.matrix(eval(parse(text=paste0("~",paste0(full_model_list, collapse="+")))), data=design_df) 

	dataObject = edgeR::DGEList(counts=count_df)
	dataObject = edgeR::calcNormFactors(dataObject, method="TMM")
	
	libsize = dataObject$samples$lib.size
	normfactors = dataObject$samples$norm.factors
	libsize = libsize * normfactors
	tmm = count_df / libsize * 1e6

	#cpm <- data.frame(edgeR::cpm(dataObject))
	#colnames(cpm) <- colnames(count_df)

	Voom = limma::voom(dataObject, design=model_matrix)

	return(list(TMM=tmm, Voom=Voom))
}

surrogate_variable_analysis <- function(normalized_count_df, design_df){

	## compute residual variance after regressing out covariables
	voom_residuals = t(resid(lm(as.formula(paste0('t(normalized_count_df)', "~", paste0(full_model_list, collapse="+"))), data=design_df)))

	### estimate N SVs
	isva_result = isva::EstDimRMT(voom_residuals, F)
	num_svs <<- isva_result$dim + 1

	### try if the model is converging with numSVs SVs, if not, rerun with 1 SV less 
	full_model_matrix = model.matrix(eval(parse(text=paste0("~",paste0(full_model_list, collapse="+")))), data=design_df) 
	null_model_matrix = model.matrix(eval(parse(text=paste0("~",config$null_model))), data=design_df)
	res = try(SmartSVA::smartsva.cpp(normalized_count_df, full_model_matrix, mod0=null_model_matrix, n.sv=num_svs, alpha=1, B=200, VERBOSE=F))
	if (inherits(res, "try-error")) {
		while (inherits(res, "try-error")){
			num_svs <<- num_svs - 1
			if(num_svs < 2) { stop("SVA model is not converging") }
			res = try(SmartSVA::smartsva.cpp(normalized_count_df, full_model_matrix, mod0=null_model_matrix, n.sv=num_svs, alpha=1, B=200, VERBOSE=F))
		}
	}
	cat(paste0(" # model is converging with estimation of ", num_svs," SVs\n"), file=log_file)
	
	## run sva !
	sv_object = SmartSVA::smartsva.cpp(normalized_count_df, full_model_matrix, mod0=null_model_matrix, n.sv=num_svs, alpha=1, B=200, VERBOSE=F)
	return(sv_object)
}

principal_component_analysis <- function(df) {

	PCA = prcomp(t(df))
	PCs = data.frame(PCA$x[, 1:num_svs])
	
	var = summary(PCA)$importance[2,1:num_svs]
	Var = data.frame("x"=1:length(var), "var"=as.vector(var)*100)
	
	if (require(ggplot2) & require(ggpubr)){
		pdf(paste0(output_path, "/estimated_variables.PCA.pdf"), width=10, height=4)
		scree = ggplot(Var, aes(x=x, y=var)) + geom_bar(stat="identity", fill="black") + ggtitle("Scree plot") + xlab("PCs") + ylab("Variance (%)")
		pc1_2 = ggplot(PCs, aes(x=PC1, y=PC2)) + geom_point(alpha=0.5) + ggtitle("PCA") + xlab(paste0("PC1 (", Var$var[1] ," %)")) + ylab(paste0("PC2 (", Var$var[2] ," %)")) + theme(legend.position="none")
		print(ggarrange(scree,pc1_2, nrow=1, widths=c(1:1)))	
		dev.off()
	}
	
	return(PCs)
}

latent_variable_evaluation <- function(voom, design, type="SV"){
	master = cbind(design, t(voom))
	
	for_correlation = cbind(design[,paste0(type, 1:num_svs)], t(voom))
	correlation = cor(for_correlation)
	correlation = correlation[-c(1:num_svs), 1:num_svs]

	if (!is.na(config$gene_annotations)){
		annot = read.csv(config$gene_annotations, sep="\t")
		rownames(annot) = annot$ID
		correlation_to_ouput = merge(annot, correlation, by="row.names")
	} else {correlation_to_ouput = correlation}
	write.table(correlation_to_ouput, file=paste0(output_path, "/estimated_variable.",type,".correlations.tsv"), sep="\t", row.names=F, quote=F)

	if (require(ggplot2) & require(reshape2) & require(ggpubr)){
		pdf(paste0(output_path, "/estimated_variables.",type,".evaluation.pdf"), width=10, height=4)
		for ( i in 1:num_svs){

			v = paste0(type, i)
			top_correlation = correlation[order(abs(correlation[,v]), decreasing=T), ][-1, ]
			top_correlation = data.frame(ID=rownames(top_correlation), gene=1:nrow(top_correlation), V=abs(top_correlation[,i]))
			top_genes = top_correlation$ID[1:3]

			subdf = master[,c(v, top_genes)]
			if(exists("annot", mode="list")){
				top_gene_names = annot[top_genes, ]$Name
				colnames(subdf) = c(v, top_gene_names)
			} else {top_gene_names = top_genes}

			melted = melt(subdf, id=c(v))
			a = ggplot(top_correlation, aes(x=gene, y=V)) + geom_point(alpha=0.5) + ylab("pearson's correlation")
			c = ggplot(melted, aes(x=melted[, v], y=value, group=variable)) + geom_point(alpha=0.5) + facet_wrap(~variable, scale="free") + theme(legend.position="bottom") + xlab(v) + ylab("voom level")
			print(ggarrange(a,c, nrow=1, widths=c(1:3)))
		}
		dev.off()
	}
}


create_corrplot = function(df, type) {
	correlation_matrix = cor(na.omit(df))
	write.table(cbind(ID=rownames(correlation_matrix),correlation_matrix), paste0(output_path, "/correlation_matrix_",type,".tsv"), sep="\t",quote=F, row.names=F)
	
	if (require(corrplot)){
		pdf(paste0(output_path, "/correlation_plot.",type,".pdf"), width=12, height=12)
		corrplot(correlation_matrix, tl.col="black",method = 'color')
		dev.off()
	}
}


run_limma <- function(count_df, design_df){

	### create covariate matrix
	model_matrix = model.matrix(eval(parse(text=paste0("~",paste0(full_model_list, collapse="+")))), data=design_df)
	### run voom with full model
	normalization_obj = edgeR_normalization(count_df, design_df)
	voom_obj = normalization_obj$Voom
	voom_df = normalization_obj$E
	#write.table(cbind("ID"=rownames(voom_df), voom_df), file=paste0(output_path, "/voom_counts.tsv" ), sep="\t", quote=F, row.names=F)

	### run the regressions
	fit <- limma::lmFit(voom_obj, model_matrix)
	contrast = colnames(fit$coefficients)[2]
	if (contrast != config$contrast) {stop(paste0("limma is using ",contrast," as contrast"))} 
	### run eBayes() to borrow information across genes
	ebayes_fit <- limma::eBayes(fit)
	save(fit, ebayes_fit, file=paste0(output_path,"/ebayes_fit.RData"))

	### Extract limma results
	results <- limma::topTable(ebayes_fit, coef=contrast, number=Inf, sort.by="p")
	results <- cbind(ID=rownames(results), results)

	### Annotate limma results
	average_tmm <- data.frame(ID=rownames(normalization_obj$TMM), avgTMM=rowMeans(normalization_obj$TMM))
	results <- merge(results, average_tmm, by="ID", all.x=T)

	if (!is.na(config$gene_annotations)) {
		gene_description = read.csv(config$gene_annotations, sep="\t")
		results <- merge(gene_description, results, by="ID", all.y=T)
	}

	results <- results[order(results$P.Value), ]
	write.table(results, file=paste0(output_path,"/limma_results.tsv"), sep="\t", quote=F, row.names=F)

}


run_deseq <- function(count_df, design_df){

	# create model for deseq where contrast is the last variable
	i = grep(config$contrast, full_model_list)
	no_contrast_model_list = full_model_list[-i]
	full_model_formula = paste0("~",paste0(no_contrast_model_list, collapse="+"),"+",config$contrast)
	
	ddsFullCountTable = DESeq2::DESeqDataSetFromMatrix(countData=count_df, colData=design_df, design=as.formula(full_model_formula))
	dds = DESeq2::DESeq(ddsFullCountTable)
	contrast = DESeq2::resultsNames(dds)[2]
	cat(paste0("DESeq2 is using : ",contrast, " as contrast !"), file=log_file)
	
	deseq_counts <- DESeq2::counts(dds, normalized=TRUE)
	deseq_counts <- cbind(ID=rownames(deseq_counts), deseq_counts)
	write.table(deseq_counts, file=paste0(output_path, "/deseq_normalized_counts.tsv" ), sep="\t", quote=F, row.names=F)
	
	save(dds, file=paste0(output_path,"/deseq2_dds_fit.RData"))
	results <- data.frame(DESeq2::results(dds, name=config$contrast))
	results <- cbind(ID=rownames(results), results)	

	if (! is.na(config$gene_annotations)) {
		gene_description = read.csv(config$gene_annotations, sep="\t")
		results <- merge(gene_description, results, by="ID", all.y=T)
		colnames(results) <- c("ID",colnames(gene_description)[-1],"baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
	} else {
		colnames(results) <- c("ID","baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val")
	}

	results <- results[order(results$P.Value),]
	write.table(results, file=paste0(output_path,"/deseq_results.tsv"), sep="\t", quote=F, row.names=F)
}


volcano_plot <- function(count_df, design_df){


}