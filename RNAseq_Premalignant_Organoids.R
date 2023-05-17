#刘虹余_RENEW_DESEQ2_premalignantIZE_premalignant
project <- c("premalignant")
deal_design <- c("TPA","TP")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- paste(getwd(),"/0_DESeq2_2vs2_premalignant_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="./mm10_anno/ebg_mm10.RData")
load(file="./mm10_anno/txdb_mm10.RData")
deal_TPA <- deal_design[1]
deal_TP <- deal_design[2]
deal_names <- paste(deal_TPA,deal_TP,sep="_VS_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"0_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"1_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"1_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"1_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"2_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"3",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"4",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"5",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"6",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"7",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"8",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"9",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})


indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "2vs2_sampletable.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",
singleEnd=FALSE,
ignore.strand=TRUE,
fragments=FALSE)
setwd(file_path)
save(se,file=se.save_RData)

colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
premalignantized_DEseq <- counts(DESeq_counts, premalignantized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","TPA","TP"))
res_1 <- cbind(premalignantized_DEseq,DEseq_res)

if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
}

res_1$ENSEMBL <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("TPA_1","TPA_2","TP_1","TP_2","DESeq2_TPA_1","DESeq2_TPA_2","DESeq2_TP_1","DESeq2_TP_2", "baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","entrez","GENENAME")
write.csv(all_summry, "premalignant_DEseq2premalignantized_TPA_VS_TP_allsummry.csv")
