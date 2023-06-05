#--------------------------------
# Differential expression analysis with SVA and limma-voom  
# by Frederique White
# 2022-11-24 
# 
# Usage : Rscript de-analysis.R -c path_to_config_file -o output_directory
#
# 
# Libraries
suppressMessages(library(rjson))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
set.seed(19) 

#--------------------------------

if(!exists("read_command_line_arguments", mode="function")) source("core.R")

## init global variables and get arguments
read_command_line_arguments()
## read config file and verify if all provided files exist
read_config_file()

obj = prepare_dataframes("dea")

clean_design_df = obj$design_df
clean_count_df = obj$count_df
# also create global variables: 
# 	- included_sample_df, 
#	- full_model_list


run_limma(clean_count_df, clean_design_df)

#run_deseq(clean_count_df, clean_design_df)
