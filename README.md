# Differential Expression analysis with limma


## Description
Differential expression analysis for RNA-seq data. From raw quantification and TPM table, the pipeline apply custom filters. Genes with high enough expression levels are used to compute surrogate variables to be included in DE analysis. DE analysis is then conducted by limma. 

## Installation
You'll need R and the following packages:
- [ ] edgeR
- [ ] isva
- [ ] SmartSVA
- [ ] limma
- [ ] ggplot2
- [ ] ggpubr
- [ ] ggrepel


## Usage
```bash
Rscript de-analysis.R -c <path_to_config_file> -o <output_directory>
```

### config file example

|variable_name|variable_type|variable_value|
|:------------|:-----------:|:-------------|
included_samples|	file|	/absolute/path/to/sample_list.txt
read_count_matrix|	file|	/absolute/path/to/gene_reads.gct.gz
tpm_count_matrix|	file|	/absolute/path/to/gene_tpm.gct.gz
design|	file|	/absolute/path/to/design.tsv
gene_annotations|	file|	/absolute/path/to/gene_annotations.tsv
mappability_scores|	file|	/absolute/path/to/gencode.v30.GRCh38.mappability.txt
clean_counts|    file|    /absolute/path/to/clean_counts.tsv
count_threshold|	parameter|	6
tpm_threshold|	parameter|	0.5
sample_frac_threshold|	parameter|	0.2
mappability_threshold|	parameter|	0.8
contrast|	parameter|	status_diabetes
model|	parameter|	status_diabetes+Sex+Age
run_SVA|	parameter|	FALSE
number_surrogate_variables|	parameter|	30
