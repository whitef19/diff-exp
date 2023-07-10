# Differential Expression analysis with limma


## Description
Differential expression analysis for RNA-seq data. From raw quantification and TPM table, the pipeline apply custom filters. Genes with high enough expression levels are used to compute surrogate variables to be included in DE analysis. DE analysis is then conducted by limma. 


## Requirements
You'll need R and the following packages:
- [ ] rjson
- [ ] edgeR
- [ ] isva
- [ ] SmartSVA
- [ ] limma

Optionally if you want to produce plots:
- [ ] ggplot2
- [ ] reshape2
- [ ] ggpubr
- [ ] ggrepel


## Usage
```bash
# run surrogate variables analysis 
Rscript sv-analysis.R -d <path_to_design_file> -c <path_to_config_file> -o <output_directory>

# run differential expression analysis 
Rscript de-analysis.R -d <path_to_design_file> -c <path_to_config_file> -o <output_directory>

# produce volcano plots 
Rscript volcano.R <results_directory>
```



## config file

The config file must contains all options listed in the **config.csv** file included in this git. Unused variables should be given the value NA.   

#### Example of config file

 |included_samples|	/absolute/path/to/included_sample.txt
 |read_count_matrix|/absolute/path/to/gene_reads.gct.gz
 |count_threshold|	6
 |normalized_count_matrix|	NA
 |tpm_threshold|	0.5
 |sample_frac_threshold|	0.2
 |mappability_scores|	NA
 |mappability_threshold|	NA
 |gene_annotations|	/absolute/path/to/gene_annotations.tsv
 |contrast|	status_diabetes
 |full_model|status_diabetes+sex+age+sequencing_run
 |null_model|sequencing_run
 |estimated_variables_method|SV
 |number_estimated_variables_to_include|0

## Precisions on files to provide

#### included_samples (required)

The file *included_samples* contains two columns with a header. The file lists the sample IDs to include in the analysis, one ID per row. The sample IDs listed in the first column (design_id) must be compatible with the IDs in *design*. The sample IDs listed in the second column (quantification_table_id) must be compatible with the IDs in *read_count_matrix* (and *normalized_count_matrix*).

```
design_id   quantification_table_id
sample_1	rnaseq1
sample_2	rnaseq2
sample_3	rnaseq3
sample_4	rnaseq4
```

#### read_count_matrix (required) and normalized_count_matrix (optional)

The file *read_count_matrix* should contains raw counts and *normalized_count_matrix* should contains normalized counts such as TPM, FPKM or RPKM. 
The rows must contain the genes and the columns must contain the samples.

The file *normalized_count_matrix* is optional. *read_count_matrix* and *normalized_count_matrix* must have the same sample IDs.

```
Name                Description     rnaseq1        rnaseq2
ENSG00000223972.5   GENE1           3               0
ENSG00000227232.5   GENE2           184             150
ENSG00000278267.1   GENE3           400             534
ENSG00000243485.5   GENE4           84              544

```

#### design (required)

The file *design* should contains all covariables used in *model*. The rows must contain the samples and the columns must contain the covariables.
(If you want to provided your own surrogate variables, add them to the design fill and name them SV1, SV2, SV3, ect.)

```
ID          status_diabetes     sex     age     SV1          SV2
sample_1    0                   0       25      0.122        0.037      
sample_2    0                   1       28      0.109        0.105  
sample_3    1                   0       23      0.092        0.007  
sample_4    1                   1       24      0.137        0.009  

```

#### mappability_scores (optional)

The file *mappability_scores* is optional. 
It should contains the mappability scores. This file should have two columns: gene IDs (compatible with the IDs in *read_count_matrix* and *normalized_count_matrix*) and mappability scores without header.

```
ENSG00000223972.5	0.626473
ENSG00000227232.5	0.152829
ENSG00000278267.1	0.372085
ENSG00000243485.5	0.908104
```

#### gene_annotations (optional)

The file *gene_annotations* is optional. It should contains a first column with the gene IDs (compatible with the IDs in *read_count_matrix* and *normalized_count_matrix*). It has a header with the first column named *ID*. You can provide any useful informations in the following columns that you want to retrieve in the final result table.

```
ID                  Name    mappability_score   annotations
ENSG00000223972.5   GENE1   0.626473            associated_with_diabetes 
ENSG00000227232.5   GENE2   0.152829            unknown  
ENSG00000278267.1   GENE3   0.372085            associated_with_insuline 
ENSG00000243485.5   GENE4   0.908104            unknown     
```
