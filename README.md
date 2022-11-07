# Differential Expression analysis with limma


## Description
Differential expression analysis for RNA-seq data. From raw quantification and TPM table, the pipeline apply custom filters. Genes with high enough expression levels are used to compute surrogate variables to be included in DE analysis. DE analysis is then conducted by limma. 


## Requirements
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

#### Example of config file

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
model|	parameter|	status_diabetes+sex+age
run_SVA|	parameter|	FALSE
number_surrogate_variables|	parameter|	2




## Precisions on files to provide

#### included_samples (required)

The file *included_samples* contains one column without header. The file lists the sample IDs to include in the analysis, one ID per row. The sample IDs listed must be compatible with the IDs in *read_count_matrix*, *normalized_count_matrix* and *design*.   

```
sample_1
sample_2
sample_3
sample_4
```

#### read_count_matrix (required) and normalized_count_matrix 

The file *read_count_matrix* should contains raw counts and *normalized_count_matrix* should contains normalized counts such as TPM, FPKM or RPKM. 
The rows must contain the genes and the columns must contain the samples.

The file *normalized_count_matrix* is optional. *read_count_matrix* and *normalized_count_matrix* must have the same sample IDs.

```
Name                Description     sample_1        sample_2
ENSG00000223972.5   GENE1           3               0
ENSG00000227232.5   GENE2           184             150
ENSG00000278267.1   GENE3           400             534
ENSG00000243485.5   GENE4           84              544

```

#### design (required)

The file *design* should contains all covariables used in *model*. The rows must contain the samples and the columns must contain the covariables.
If you want to provided your own surrogate variables, add them to the design fill and name them SV1, SV2, SV3, ect. 

```
ID          status_diabetes     sex     age     SV1          SV2
sample_1    0                   0       25      0.122        0.037      
sample_2    0                   1       28      0.109        0.105  
sample_3    1                   0       23      0.092        0.007  
sample_4    1                   1       24      0.137        0.009  

```

#### mappability_scores

The file *mappability_scores* is optional. 
It should contains the mappability scores. This file should have two columns: gene IDs (compatible with the IDs in *read_count_matrix* and *normalized_count_matrix*) and mappability scores without a header.

```
ENSG00000223972.5	0.626473
ENSG00000227232.5	0.152829
ENSG00000278267.1	0.372085
ENSG00000243485.5	0.908104
```

#### gene_annotations

The file *gene_annotations* is optional. It should contains a first column with the gene IDs (compatible with the IDs in *read_count_matrix* and *normalized_count_matrix*). It has a header with the first column named *ID*. You can provide any useful informations in the following columns that you want to retrieve in the final result table.

```
ID                  Name    mappability_score   annotations
ENSG00000223972.5   GENE1   0.626473            associated_with_diabetes 
ENSG00000227232.5   GENE2   0.152829            unknown  
ENSG00000278267.1   GENE3   0.372085            associated_with_insuline 
ENSG00000243485.5   GENE4   0.908104            unknown     
```
