

suppressMessages(library(ggplot2))
suppressMessages(require(ggrepel))
suppressMessages(library(ggpubr))


ARG = commandArgs(trailingOnly = T)
indir = ARG[1]


results <- read.csv(paste0(indir,"/limma_results.tsv"), sep="\t", row.names=1)

config <- read.csv(paste0(indir,"/config_file_used"), sep="\t", row.names=1)
model <- config["full_model","variable_value"]
nbsv <- config["number_surrogate_variables","variable_value"]

design_file <- config["design","variable_value"]
design <- read.csv(design_file, sep="\t", row.names=1)

nbsvs = colnames(design)[ncol(design)]
nbsamples = nrow(design)
results$absFC = abs(results$logFC)
logFC_threshold <- (mean(results$absFC) + 5 * sd(results$absFC))

top_hits <- subset(results, absFC > logFC_threshold)
top_hits <- subset(results, P.Value < 0.05 )
results$color <- ifelse(results$absFC > logFC_threshold, "red", "grey")
results$color <- ifelse(results$P.Value > 0.001, "grey", results$color)


pdf(paste0(indir, "/redone.volcano.limma.pdf"), width=7, height=6)

volcano <- ggplot(results, aes(x=logFC, y=-log10(P.Value), color=color)) + geom_point(alpha=0.5) + 
	ylab("-log10 nominal p-value") + xlab("log2 FC") + 
	xlim(-max(results$absFC), max(results$absFC)) +
	ggtitle(paste0("~ ", model, "\n + ", nbsvs,"\n n = ", nbsamples, " samples")) + 
	geom_vline(xintercept=c(-logFC_threshold, logFC_threshold), color="grey", linetype="dashed") +
	geom_hline(yintercept=c(-log10(0.05), -log10(0.001)), color="grey",linetype="dashed") +
	#geom_hline(yintercept=-log10(0.05), color="grey", linetype="dashed") +
	theme(plot.title=element_text(size=5), legend.position="none") + scale_color_identity()

if (nrow(subset(results, absFC > logFC_threshold & P.Value < 0.05)) > 0) {
	if ("Name" %in% names(results)) {
		print(volcano + geom_text_repel(data =subset(results, P.Value < 0.001), size=3, aes(label=Name), show.legend=F, colour = "black"))
	} else {
		print(volcano + geom_text_repel(data=subset(results,  P.Value < 0.001), aes(label=ID), show.legend=F, colour = "black"))
	}
} else {print(volcano)}

dev.off()