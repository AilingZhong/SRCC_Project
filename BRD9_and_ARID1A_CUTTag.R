************************************Arid1a结合区域分析***********************************************************
************************************Arid1a结合区域分析***********************************************************
************************************Arid1a结合区域分析***********************************************************
************************************Arid1a结合区域分析***********************************************************
#q=0.00001

macs2 callpeak -t TP1_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TP1_Brd9_CUTTag

macs2 callpeak -t TP2_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TP2_Brd9_CUTTag

macs2 callpeak -t TP3_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TP3_Brd9_CUTTag


macs2 callpeak -t TPA1_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TPA1_Brd9_CUTTag

macs2 callpeak -t TPA2_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TPA2_Brd9_CUTTag

macs2 callpeak -t TPA3_Brd9_CUTTag.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Brd9_CUTTag_q_e05 \
-n TPA3_Brd9_CUTTag


macs2 callpeak -t Arid1a_CUTTAG_TP2.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Arid1a_CUTTag_q_e05 \
-n anti_Arid1a_2_CUTTag

macs2 callpeak -t Arid1a_CUTTAG_TP1.filter_dupli_chrM_last.bam \
-q 0.00001 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir Arid1a_CUTTag_q_e05 \
-n anti_Arid1a_1_CUTTag


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicFeatures)
library(ReactomePA)
library(AnnotationDbi)
library(DOSE)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

TP1_Arid1a <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TP1_Arid1a.bed")
TP2_Arid1a <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TP2_Arid1a.bed")
TP1_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TP1_Brd9.bed")
TP2_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TP2_Brd9.bed")
TP3_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TP3_Brd9.bed")
TPA1_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TPA1_Brd9.bed")
TPA2_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TPA2_Brd9.bed")
TPA3_Brd9 <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/1_LHY/4_Brd9/fastp_QC/filter_bam/filter_MT_bam/TPA_VS_TP_BRD9/TPA3_Brd9.bed")


**************************************试试用macs2定量比较差异peak*********************************************************
**************************************试试用macs2定量比较差异peak*********************************************************
**************************************试试用macs2定量比较差异peak*********************************************************
**************************************试试用macs2定量比较差异peak*********************************************************
**************************************试试用macs2定量比较差异peak*********************************************************
**************************************试试用macs2定量比较差异peak*********************************************************
XY_runConsensusRegions <- function (testRanges, method = "majority", overlap = "any")
{
    if (length(testRanges) > 1) {
        reduced <- reduce(unlist(testRanges))
        consensusIDs <- paste0("consensus_", 
seq(1, length(reduced)))
        mcols(reduced) <- do.call(cbind, lapply(testRanges, function(x) (reduced %over%
            x) + 0))
        if (method == "majority") {
            reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) >
                length(testRanges)/2, ]
        }
        if (method == "none") {
            reducedConsensus <- reduced
        }
        if (is.numeric(method)) {
            reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) >
                method, ]
        }
        consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
        mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)),
            consensusIDs)
        return(reducedConsensus)
    }
}
environment(XY_runConsensusRegions) <- asNamespace('soGGi')


#这是抽取他们merge的所有common的peaks开始定量计算每一个region的counts
library(dplyr)
library(limma)
library(trqwe)
files <- list.files(path = "./Brd9_CUTTag_q_e05", pattern = "peaks.narrowPeak$", full.names = TRUE)
name <- gsub("_Brd9_CUTTag_peaks.narrowPeak","",basename(files))
all_myPeaks <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
names(all_myPeaks) <- name
myPeaks <- all_myPeaks[c("TPA1","TPA2","TPA3","TP1","TP2","TP3")]

consensusToCount <- XY_runConsensusRegions(GRangesList(myPeaks), "none")
consensusToCount
mcsaveRDS(consensusToCount, file = "Brd9_binding_TPA_VS_TP_consensusToCount.rds",mc.cores=20)
library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% rowSums
table(occurrences) %>% rev %>% cumsum
bamsToCount <- dir("./filter_MT_bam/TPA_VS_TP_BRD9", full.names = TRUE, pattern = "*.\\.filter_dupli_chrM_last.bam$")
# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
    start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
    Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
    countMultiMappingReads = FALSE, maxFragLength = 100)
myCounts <- fcResults$counts
colnames(myCounts) <- names(all_myPeaks)
mcsaveRDS(myCounts, file = "Brd9_binding_TPA_VS_TP_countsFromBAM.rds",mc.cores=20)


consensusToCount <- readRDS("./Brd9_binding_TPA_VS_TP_consensusToCount.rds")
myCounts <- readRDS("./Brd9_binding_TPA_VS_TP_countsFromBAM.rds")

sel_myCounts <- myCounts[,c("TPA1","TPA2","TPA3","TP1","TP2","TP3")]
library(DESeq2)
Group <- factor(c("TPA","TPA","TPA","TP","TP","TP"))
metaData <- data.frame(Group, row.names = colnames(sel_myCounts))
metaData$type <- rownames(metaData)
library(DESeq2)
TPA_VS_TP_Brd9_DDS <- DESeqDataSetFromMatrix(sel_myCounts, metaData, ~Group, rowRanges = consensusToCount)
TPA_VS_TP_Brd9_DDS <- DESeq(TPA_VS_TP_Brd9_DDS,parallel=TRUE)
TPA_VS_TP_Brd9_DDS_sel <- TPA_VS_TP_Brd9_DDS[,c("TPA1","TPA2","TPA3","TP1","TP2","TP3")]
pri_data <- counts(TPA_VS_TP_Brd9_DDS_sel, normalized=TRUE)
pri_data <- data.frame(pri_data)
pri_data$region <- rownames(pri_data)
library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)
TPA_VS_TP_Brd9_DDS_sel <- TPA_VS_TP_Brd9_DDS[,c("TPA1","TPA2","TPA3","TP1","TP2","TP3")]
TPA_VS_TP_Brd9_data <- counts(TPA_VS_TP_Brd9_DDS_sel, normalized=TRUE)
TPA_VS_TP_Brd9_data <- data.frame(TPA_VS_TP_Brd9_data)
TPA_VS_TP_Brd9_data$region <- rownames(TPA_VS_TP_Brd9_data)
TPA_VS_TP_Brd9_DDS_results <- results(TPA_VS_TP_Brd9_DDS, c("Group", "TPA","TP"), format = "GRanges",parallel=TRUE)
TPA_VS_TP_Brd9_DDS_results <- TPA_VS_TP_Brd9_DDS_results[order(TPA_VS_TP_Brd9_DDS_results$pvalue)]
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno_TPA_VS_TP_Brd9_DDS_results <- annotatePeak(TPA_VS_TP_Brd9_DDS_results, tssRegion=c(-3000, 3000),TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
anno_TPA_VS_TP_Brd9_DDS_results <- data.frame(anno_TPA_VS_TP_Brd9_DDS_results)
anno_TPA_VS_TP_Brd9_DDS_results$region <- paste(anno_TPA_VS_TP_Brd9_DDS_results$seqnames,sep="_",anno_TPA_VS_TP_Brd9_DDS_results$start)
anno_TPA_VS_TP_Brd9_DDS_results$region <- paste("ID",sep="_",anno_TPA_VS_TP_Brd9_DDS_results$region)
anno_TPA_VS_TP_Brd9_DDS_results$region <- paste(anno_TPA_VS_TP_Brd9_DDS_results$region,sep="_",anno_TPA_VS_TP_Brd9_DDS_results$end)

all_file <- merge(pri_data,anno_TPA_VS_TP_Brd9_DDS_results,by="region")
write.csv(all_file,"1_TPA_VS_TP_Brd9_all_anno_peaks.csv")
p005_peaks <- subset(all_file,pvalue<0.05)
write.csv(p005_peaks,"2_TPA_VS_TP_Brd9_peaks_p005.csv")