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

## config file

The config file must contains all options listed in the **config.tsv** file included in this git. Unused variables should be given the value NA.   


### Example

|variable_name|variable_type|variable_value|
|:------------|:-----------:|:-------------|
included_samples|	file|	/absolute/path/to/sample_list.txt
read_count_matrix|	file|	/absolute/path/to/gene_reads.gct.gz
normalized_count_matrix|	file|	NA
design|	file|	/absolute/path/to/design.tsv
gene_annotations|	file|	/absolute/path/to/gene_annotations.tsv
mappability_scores|	file|	NA
clean_counts|    file|    /absolute/path/to/clean_counts.tsv
count_threshold|	parameter|	6
tpm_threshold|	parameter|	0.5
sample_frac_threshold|	parameter|	0.2
mappability_threshold|	parameter|	NA
contrast|	parameter|	status_diabetes
model|	parameter|	status_diabetes+Sex+Age
run_SVA|	parameter|	FALSE
number_surrogate_variables|	parameter|	30

## Precisions on files to provide

### included_samples

The file *included_samples* contains one column without header. The file lists the sample IDs to include in the analysis, one ID per row. The sample IDs listed must be compatible with the IDs in *read_count_matrix*, *normalized_count_matrix* and *design*.   

```
sample_1
sample_2
sample_3
sample_4
```

### read_count_matrix and normalized_count_matrix

The file *read_count_matrix* should contains raw counts and *normalized_count_matrix* should contains normalized counts such as TPM, FPKM or RPKM. 
Rows contain genes and columns contain samples.  

```
Name	            Description	sample_1        sample_2
ENSG00000223972.5	GENE1 	    3               0
ENSG00000227232.5	GENE2       184             150
ENSG00000278267.1	GENE3	    400             534
ENSG00000243485.5	GENE4	    84              544

```

### design


```
ID          status_diabetes	sex	age
sample_1	0	            0	25
sample_2	0	            1	28
sample_3	1	            0	23
sample_4	1	            1	24

```

### mappability_scores


```
ENSG00000223972.5	0.626473
ENSG00000227232.5	0.152829
ENSG00000278267.1	0.372085
ENSG00000243485.5	0.908104
```
