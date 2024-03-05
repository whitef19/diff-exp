#--------------------------------
# Surrogate variables analysis for differential expression analysis  
# by Frederique White
# 2023-05-25 
# 
# Usage : Rscript sv-analysis.R -c path_to_config_file -d path_to_design_file -o output_directory
# 
# Libraries

suppressMessages(library(rjson))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(isva))
suppressMessages(library(SmartSVA))

set.seed(48121620) 

#--------------------------------

if(!exists("read_command_line_arguments", mode="function")) source("core.R")

## init global variables and get arguments from command line
read_command_line_arguments()

## read config file and verify if all provided files exist
read_config_file()

### read "included_samples"

obj = prepare_sva_dataframes()
clean_design_df = obj$design_df
clean_count_df = obj$count_df
# also create global variables: 
# 	- included_sample_df, 
#	- full_model_list

## write clean dataframe:
write.table(cbind("ID"=rownames(clean_count_df), clean_count_df), file=paste0(output_path,"/clean_gene_counts.tsv"), sep="\t", quote=F, row.names=F)



### Surrogate Variable Analysis

## compute voom normalized counts for SVA
normalization_obj = edgeR_normalization(clean_count_df, clean_design_df)

sv_object = surrogate_variable_analysis(normalization_obj$Voom$E, clean_design_df)
# also create global variables:
# 	- num_svs
SVs = as.data.frame(sv_object$sv)
names(SVs) = paste0("SV",1:ncol(SVs))
## concatenate SVs to covariate dataframe 
clean_design_df = cbind(clean_design_df, SVs)

# if PCA is more appropriate :
### Principal Component Analysis
#PCs = principal_component_analysis(normalization_obj$Voom$E)
## concatenate PCs to covariate dataframe 
#clean_design_df = cbind(clean_design_df, PCs)


write.table(cbind("ID"=rownames(clean_design_df), clean_design_df), file=paste0(output_path,"/full_design.tsv"), sep="\t", row.names=F, quote=F)

### evaluate SVs 
included_covariate_list = c(full_model_list, names(SVs))
create_corrplot(clean_design_df[, included_covariate_list], "SV")
latent_variable_evaluation(normalization_obj$Voom$E, clean_design_df[, included_covariate_list], "SV")

#included_covariate_list = c(full_model_list, names(PCs))
#create_corrplot(clean_design_df[, included_covariate_list], "PC")
#latent_variable_evaluation(normalization_obj$Voom$E, clean_design_df[, included_covariate_list], "PC")

#included_covariate_list = c(full_model_list, names(SVs), names(PCs))
#create_corrplot(clean_design_df[, included_covariate_list], "SV_and_PC")

cat(paste0("##### All done !: ",Sys.Date(), " at ", Sys.time(),"\n"), file=log_file)
