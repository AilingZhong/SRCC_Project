*******************************dBRD9_VS_DMSO_RNAseq*********************************************************
*******************************dBRD9_VS_DMSO_RNAseq*********************************************************
*******************************dBRD9_VS_DMSO_RNAseq*********************************************************
*******************************dBRD9_VS_DMSO_RNAseq*********************************************************
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

project <- c("SRCC")
deal_design <- c("inhibitor","DMSO")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- paste(getwd(),"/3_inhibitor_VS_DMSO_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="./mm10_anno/ebg_mm10.RData")
load(file="./mm10_anno/txdb_mm10.RData")
deal_inhibitor <- deal_design[1]
deal_DMSO <- deal_design[2]
deal_names <- paste(deal_inhibitor,deal_DMSO,sep="_VS_")
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
csvfile <- file.path(indir, "sampleTables.csv")
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

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL

pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()


pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()


colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL



pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()


pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()


colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","inhibitor","DMSO"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

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
res_1$ENTREZID <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("inhibitor_1","inhibitor_2","inhibitor_3","DMSO_1","DMSO_2","DMSO_3","DESeq2_inhibitor_1","DESeq2_inhibitor_2","DESeq2_inhibitor_3","DESeq2_DMSO_1","DESeq2_DMSO_2", "DESeq2_DMSO_3","baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "SYMBOL","ENTREZID","GENENAME")
write.csv(all_summry, "renew_DEseq2_inhibitor_VS_DMSO_allsummry.csv")