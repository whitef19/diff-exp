

library(ggplot2)
library(ggpubr)
library(reshape2)

#indir = "outputs/maternal/oldfashion/sva/"
#nbsv = 37
ARG = commandArgs(trailingOnly = T)
indir = ARG[1]
nbsv = ARG[2]


df <- read.csv(paste0(indir, "/voom_counts.tsv"), sep="\t", row.names=1, check.names=F)
data <- as.data.frame(t(df))

design <- read.csv(paste0(indir, "/design_svs.tsv"),sep="\t", row.names=1)
df <- merge(design, data, by="row.names")

annot <- read.csv("data/gene_annotations.tsv", sep="\t")
rownames(annot) <- annot$ID

for_correlation <- merge(design[,paste0("SV", 1:nbsv)], data, by="row.names")
correlation <- cor(for_correlation[, -1])
correlation <- correlation[, 1:nbsv]

correlation_to_ouput <- merge(annot, correlation, by="row.names")
write.table(correlation_to_ouput, file=paste0(indir, "/sv_correlations.tsv"), sep="\t", row.names=F, quote=F)


pdf(paste0(indir, "/sva_evaluations.pdf"))
for ( i in 1:nbsv){
	SV = paste0("SV", i)
	tmp = correlation[order(abs(correlation[,SV]), decreasing=T), ][-1, ]
	svcor = data.frame(ID=rownames(tmp), gene=1:nrow(tmp), SV=abs(tmp[,i]))
	a = ggplot(svcor, aes(x=gene, y=SV)) + geom_point(alpha=0.5) + ylab("pearson correlation")
	b = ggplot(svcor, aes(x=SV)) + geom_histogram(alpha=0.5, bins=20) + xlab("pearson correlation")
	r1 = ggarrange(a,b, nrow=1)

	top_genes <- svcor$ID
	top_genes <- top_genes[!"SV" %in% top_genes]
	top_genes <- top_genes[1:2]
	top_gene_names <- annot[top_genes, ]$Name
	subdf <- df[,c("Patient_No_etude", SV,"Accouchement_Sexe", top_genes)]
	colnames(subdf) <- c("Patient_No_etude", SV, "Sex", top_gene_names)
	melted <- melt(subdf, id=c("Patient_No_etude", SV,"Sex"))
	melted$Sex <- as.factor(melted$Sex) 
	c = ggplot(melted, aes(x=melted[, SV], y=value, group=variable, color=Sex)) + geom_point(alpha=0.5) + facet_wrap(~variable, scale="free") + theme(legend.position="bottom") + xlab(SV) + ylab("voom level")
	
	svcor <- cbind(svcor[1:15,], Name=annot[svcor[1:15,]$ID, ]$Name)
	d = ggtexttable(svcor[, c("Name", "SV")], rows = NULL, theme = ttheme("light", base_size=9, padding=unit(c(2, 2),"mm")))
	r2 = ggarrange(c,d, nrow=1, widths=c(3,2))

	print(ggarrange(r1,r2, nrow=2))
}

dev.off()
