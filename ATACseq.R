****************************************************************************
****************************************************************************
******************************ATACseq_analysis******************************
****************************************************************************
****************************************************************************
****************************************************************************
vim Arid1a_samplelist
Arid1a_KO_1 ./A1_ATAC_FKDL202618100-1a
Arid1a_KO_2 ./A2_ATAC_FKDL202618101-1a
Control_1 ./S1_ATAC_FKDL202618098-1a
Control_2 ./S2_ATAC_FKDL202618099-1a

cat Arid1a_samplelist | while read id ; do
arr=($id)
fq2=${arr[1]}'_2.fq.gz'
fq1=${arr[1]}'_1.fq.gz'
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
NGmerge -a -1 $fq1 -2 $fq2 -o $sample -u 41 -v -n 30


bowtie2 -p 20 --very-sensitive -X2000 -x ./bowtie2_mm10_index/mm10 \
-1 Arid1a_KO_1_1.fastq.gz -2 Arid1a_KO_1_2.fastq.gz | samtools sort -O bam -@ 20 -o - > Arid1a_KO_1.NGmerge.sorted.bam
bowtie2 -p 20 --very-sensitive -X2000 -x ./bowtie2_mm10_index/mm10 \
-1 Arid1a_KO_2_1.fastq.gz -2 Arid1a_KO_2_2.fastq.gz | samtools sort -O bam -@ 20 -o - > Arid1a_KO_2.NGmerge.sorted.bam
bowtie2 -p 20 --very-sensitive -X2000 -x ./bowtie2_mm10_index/mm10 \
-1 Control_1_1.fastq.gz -2 Control_1_2.fastq.gz | samtools sort -O bam -@ 20 -o - > Control_1.NGmerge.sorted.bam
bowtie2 -p 20 --very-sensitive -X2000 -x ./bowtie2_mm10_index/mm10 \
-1 Control_2_1.fastq.gz -2 Control_2_2.fastq.gz | samtools sort -O bam -@ 20 -o - > Control_2.NGmerge.sorted.bam




mkdir -p filter_bam/filter_MT_bam
mkdir bed
mkdir peak_file


#####消除PCR重复######
./programme/gatk-4.1.3.0/gatk --java-options "-Xmx30G -Djava.io.tmpdir=./" MarkDuplicates \
-I Arid1a_KO_1.NGmerge.sorted.bam \
-O ./filter_bam/Arid1a_KO_1.filtdup.bam \
-M ./filter_bam/Arid1a_KO_1.dups.txt \
-REMOVE_DUPLICATES=true

./programme/gatk-4.1.3.0/gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" MarkDuplicates \
-I Arid1a_KO_2.NGmerge.sorted.bam \
-O ./filter_bam/Arid1a_KO_2.filtdup.bam \
-M ./filter_bam/Arid1a_KO_2.dups.txt \
-REMOVE_DUPLICATES=true

./programme/gatk-4.1.3.0/gatk --java-options "-Xmx30G -Djava.io.tmpdir=./" MarkDuplicates \
-I Control_1.NGmerge.sorted.bam \
-O ./filter_bam/Control_1.filtdup.bam \
-M ./filter_bam/Control_1.dups.txt \
-REMOVE_DUPLICATES=true

./programme/gatk-4.1.3.0/gatk --java-options "-Xmx30G -Djava.io.tmpdir=./" MarkDuplicates \
-I Control_2.NGmerge.sorted.bam \
-O ./filter_bam/Control_2.filtdup.bam \
-M ./filter_bam/Control_2.dups.txt \
-REMOVE_DUPLICATES=true



cat Arid1a_samplelist | while read id ; do
arr=($id)
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
samtools index  ./filter_bam/$sample.filtdup.bam
mtReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

samtools flagstat  ./filter_bam/$sample.filtdup.bam > $sample.rmdup.stat
samtools view  -h  -f 2 -q 30 ./filter_bam/$sample.filtdup.bam |grep -v chrM > ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam
samtools sort ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam  -@ 8 -o - > ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam
rm -r ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam

samtools index ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam 
samtools flagstat  ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam > $sample.filter_dupli_chrM_last.stat 
bedtools bamtobed -i ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam  > ./bed/$sample.bed ;
done 




ls *filter_dupli_chrM_last.bam |while read id;do
bamCoverage -p 20 -bs=1 --normalizeUsing BPM -b $id -o ./${id%%.*}.bs1.last.bw ;
done

for id in *.filter_dupli_chrM_last.bam; do
echo $id 
sample=$(basename $id .filter_dupli_chrM_last.bam)
echo $sample
samtools view -@ 20 -H $id | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > $sample.genome.info;
done


cat Arid1a_samplelist | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
HMMRATAC=/mnt/data/userdata/xiangyu/programme/HMMRATAC-1.2.5/HMMRATAC_V1.2.5_exe.jar
java -Xmx30g -jar $HMMRATAC -b $sample.filter_dupli_chrM_last.bam -i $sample.filter_dupli_chrM_last.bam.bai -g $sample.genome.info -o $sample ; 
done



for id in *.sorted.bam; do
echo $id 
sample=$(basename $id.sorted.bam)
echo $sample
samtools view -@ 20 -H $id | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > $sample.genome.info;
done


HMMRATAC=./programme/HMMRATAC-1.2.5/HMMRATAC_V1.2.5_exe.jar
java -Xmx30g -jar $HMMRATAC \
-b merge_Arid1a_WT.sorted.bam \
-i merge_Arid1a_WT.sorted.bam.bai \
-g merge_Arid1a_WT.sorted.bam.sorted.bam.genome.info -o merge_Arid1a_WT_HMMRATAC

HMMRATAC=./programme/HMMRATAC-1.2.5/HMMRATAC_V1.2.5_exe.jar
java -Xmx30g -jar $HMMRATAC \
-b merge_Arid1a_KO.sorted.bam \
-i merge_Arid1a_KO.sorted.bam.bai \
-g merge_Arid1a_KO.sorted.bam.sorted.bam.genome.info -o merge_Arid1a_KO_HMMRATAC



for id in *filter_dupli_chrM_last.bam;do
sample_name=${id%%.*}
echo $sample_name.length.txt
samtools view -@ 20 $id |awk '{print $9}'  > $sample_name.length.txt
done


library("future.apply")
library("ggplot2")
files <- list.files(path = "./filter_bam/filter_MT_bam", pattern = ".length.txt$", full.names = TRUE)
name <- gsub(".length.txt","",basename(files))
name <- gsub("Mmu_","",name)
all_peaks <- future_lapply(1:length(files),function(i){
  aa <- as.data.frame(data.table::fread(files[[i]]))
  colnames(aa) <- "Insertion_Size_distribution"
  aa$Insertion_Size_distribution <- abs(aa$Insertion_Size_distribution)
  aa$group <- name[i]
  return(aa)
  })
all_peaks_files <- as.data.frame(data.table::rbindlist(all_peaks))


ff1 <- ggplot(data = all_peaks_files, mapping = aes(x = Insertion_Size_distribution, fill = group)) +
            geom_histogram(bins = 200) + facet_wrap(~group, scales = "free_y") +
            xlim(c(0, 800)) + theme_bw() + theme(legend.position = "none")
ggsave("all_peaks_files.svg", plot=ff1,width = 12,height = 5,dpi=1080)





bed=/mnt/data/public_data/reference/Chip_ATACseq/mm10_bed/mm10_RefSeq.bed
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
-R $bed \
-S ATAC_merge_Arid1a_WT.bs1.bw ATAC_merge_Arid1a_KO.bs1.bw \
--numberOfProcessors 30 --skipZeros -o merge_ATAC.mat.gz 

plotHeatmap -m merge_ATAC.mat.gz -out renew_merge_ATAC.png \
--colorList  'white,red' 'white,red' \
--samplesLabel "TPS_ATAC" "TPA_ATAC" \
--plotFileFormat png



sudo mv /mnt/data/sequencedata/chip_seq/chip_seq_10_LHY_20201220_H3K27ac_4samples/rawdata /mnt/data/sequencedata/CUTTAG/1_LHY_H3k27ac_4samples
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/1.rawdata/LHY_TPA2_H3K27ac_FKDL210000783-1a/LHY_TPA2_H3K27ac_FKDL210000783-1a_2.fq.gz /mnt/data/sequencedata/CUTTAG/1_LHY_H3k27ac_4samples/rawdata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/1.rawdata/LHY_TPA2_H3K27ac_FKDL210000783-1a/LHY_TPA2_H3K27ac_FKDL210000783-1a_1.fq.gz /mnt/data/sequencedata/CUTTAG/1_LHY_H3k27ac_4samples/rawdata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/1.rawdata/LHY_TPS2_Arid1a_FKDL210000785-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/rawdata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/1.rawdata/LHY_TPS1_Arid1a_FKDL210000784-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/rawdata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/1.rawdata/IgG_FKDL210000786-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/rawdata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/2.cleandata/LHY_TPS2_Arid1a_FKDL210000785-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/cleandata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/2.cleandata/LHY_TPS1_Arid1a_FKDL210000784-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/cleandata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/2.cleandata/IgG_FKDL210000786-1a /mnt/data/sequencedata/CUTTAG/2_LHY_Arid1a_3samples/cleandata
sudo mv /mnt/data/sequencedata/CUTTAG/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J004/2.cleandata/LHY_TPA2_H3K27ac_FKDL210000783-1a /mnt/data/sequencedata/CUTTAG/1_LHY_H3k27ac_4samples/cleandata




bed=/mnt/data/public_data/reference/Chip_ATACseq/mm10_bed/mm10_RefSeq.bed
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
-R $bed \
-S Control_1.bs1.last.bw Control_2.bs1.last.bw Arid1a_KO_1.bs1.last.bw Arid1a_KO_2.bs1.last.bw \
--numberOfProcessors 30 --skipZeros -o global_matrix_TSS.mat.gz 


plotHeatmap -m global_matrix_TSS.mat.gz -out global_matrix_TSS.png \
--colorList  'white,red' 'white,red' 'white,red' 'white,red' \
--samplesLabel "Control_1" "Control_2" "Arid1a_KO_1" "Arid1a_KO_2" \
--plotFileFormat png







*************************************定量比较差异peak*****************************************************
*************************************定量比较差异peak*****************************************************
*************************************定量比较差异peak*****************************************************
*************************************定量比较差异peak*****************************************************

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
files <- list.files(path = "./filter_bam/filter_MT_bam/HMMRATAC", pattern = "_peaks.gappedPeak$", full.names = TRUE)
name <- gsub("_peaks.gappedPeak","",basename(files))
all_myPeaks <- lapply(files, ChIPQC:::GetGRanges, simple = TRUE)
names(all_myPeaks) <- name
myPeaks <- all_myPeaks[c("Arid1a_KO_1","Arid1a_KO_2","Control_1","Control_2")]
consensusToCount <- XY_runConsensusRegions(GRangesList(myPeaks), "none")
consensusToCount
mcsaveRDS(consensusToCount, file = "Arid1a_consensusToCount.rds",mc.cores=20)




library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% rowSums
table(occurrences) %>% rev %>% cumsum
bamsToCount <- dir("/mnt/data/userdata/abao/project/5_ATAC/1_LHY/NGmerge_produce_fq/filter_bam/filter_MT_bam", full.names = TRUE, pattern = "*.\\.filter_dupli_chrM_last.bam$")

# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
    start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
    Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
    countMultiMappingReads = FALSE, maxFragLength = 100)
myCounts <- fcResults$counts
colnames(myCounts) <- names(all_myPeaks)
mcsaveRDS(myCounts, file = "Arid1a_countsFromATAC.rds",mc.cores=20)





********************************Arid1a_ko_vs_control*******************************************************
********************************Arid1a_ko_vs_control*******************************************************
********************************Arid1a_ko_vs_control*******************************************************
********************************Arid1a_ko_vs_control*******************************************************
********************************Arid1a_ko_vs_control*******************************************************
********************************Arid1a_ko_vs_control*******************************************************
sel_myCounts <- myCounts[,c("Arid1a_KO_1","Arid1a_KO_2","Control_1","Control_2")]

library(DESeq2)
Group <- factor(c("Arid1a_KO","Arid1a_KO","Control","Control"))
metaData <- data.frame(Group, row.names = colnames(sel_myCounts))
metaData$type <- rownames(metaData)
library(DESeq2)
atacDDS <- DESeqDataSetFromMatrix(sel_myCounts, metaData, ~Group, rowRanges = consensusToCount)
atacDDS <- DESeq(atacDDS,parallel=TRUE)

library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)

atacDDS_sel <- atacDDS[,c("Arid1a_KO_1","Arid1a_KO_2","Control_1","Control_2")]
pri_data <- counts(atacDDS_sel, normalized=TRUE)
pri_data <- data.frame(pri_data)
pri_data$region <- rownames(pri_data)

atacDDS_results <- results(atacDDS, c("Group", "Arid1a_KO","Control"), format = "GRanges",parallel=TRUE)
atacDDS_results <- atacDDS_results[order(atacDDS_results$pvalue)]


library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
anno_atacDDS_results <- annotatePeak(atacDDS_results, tssRegion=c(-3000, 3000),TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
anno_atacDDS_results <- data.frame(anno_atacDDS_results)
anno_atacDDS_results$region <- paste(anno_atacDDS_results$seqnames,sep="_",anno_atacDDS_results$start)
anno_atacDDS_results$region <- paste("ID",sep="_",anno_atacDDS_results$region)
anno_atacDDS_results$region <- paste(anno_atacDDS_results$region,sep="_",anno_atacDDS_results$end)

all_file <- merge(pri_data,anno_atacDDS_results,by="region")
write.csv(all_file,"1_Arid1a_ATAC_all_anno_peaks.csv")

p005_peaks <- subset(all_file,pvalue<0.05)
write.csv(p005_peaks,"2_Arid1a_ATAC_peaks_p005.csv")



**********************************merge**igv**************************************************
**********************************merge**igv**************************************************
**********************************merge**igv**************************************************
**********************************merge**igv**************************************************
**********************************merge**igv**************************************************
**********************************merge**igv**************************************************

samtools merge -@ 20 Merge_ATAC_Arid1a_WT.bam Control_1.filter_dupli_chrM_last.bam Control_2.filter_dupli_chrM_last.bam
samtools merge -@ 20 Merge_ATAC_Arid1a_KO.bam Arid1a_KO_1.filter_dupli_chrM_last.bam Arid1a_KO_2.filter_dupli_chrM_last.bam

samtools sort -@ 20 Merge_ATAC_Arid1a_WT.bam -o Merge_ATAC_Arid1a_WT.sorted.bam
samtools sort -@ 20 Merge_ATAC_Arid1a_KO.bam -o Merge_ATAC_Arid1a_KO.sorted.bam

samtools index Merge_ATAC_Arid1a_WT.sorted.bam
samtools index Merge_ATAC_Arid1a_KO.sorted.bam

bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_ATAC_Arid1a_WT.sorted.bam -o Merge_ATAC_Arid1a_WT.bs1.bw
bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_ATAC_Arid1a_KO.sorted.bam -o Merge_ATAC_Arid1a_KO.bs1.bw


bed=/mnt/data/public_data/reference/Chip_ATACseq/mm10_bed/mm10_RefSeq.bed
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
-R $bed \
-S Merge_ATAC_Arid1a_WT.bs1.bw Merge_ATAC_Arid1a_KO.bs1.bw \
--numberOfProcessors 30 --skipZeros -o renew_Merge_ATAC.mat.gz 

plotHeatmap -m renew_Merge_ATAC.mat.gz -out renew_merge_ATAC.png \
--colorList "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" \
--samplesLabel "TPS_ATAC" "TPA_ATAC" \
--plotFileFormat png


computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
-R Distal_intergenic_all_peaks.bed \
-S Merge_ATAC_Arid1a_WT.bs1.bw Merge_ATAC_Arid1a_KO.bs1.bw \
--numberOfProcessors 30 --skipZeros -o distal_Merge_ATAC.mat.gz 

plotHeatmap -m distal_Merge_ATAC.mat.gz -out distal_merge_ATAC.png \
--colorList "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" \
--samplesLabel "TPS_ATAC" "TPA_ATAC" \
--plotFileFormat png


computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
-R promoter_ATAC_all_peaks.bed \
-S Merge_ATAC_Arid1a_WT.bs1.bw Merge_ATAC_Arid1a_KO.bs1.bw \
--numberOfProcessors 30 --skipZeros -o promoter_Merge_ATAC.mat.gz 

plotHeatmap -m promoter_Merge_ATAC.mat.gz -out promoter_merge_ATAC.png \
--colorList "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" "#fff5f0,#fee0d2,#fcbba1,#fc9272,#fb6a4a,#ef3b2c,#cb181d,#a50f15,#67000d" \
--samplesLabel "TPS_ATAC" "TPA_ATAC" \
--plotFileFormat png





for id in *.sorted.bam; do
echo $id 
sample=$(basename $id.sorted.bam)
echo $sample
samtools view -@ 20 -H $id | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > $sample.genome.info;
done

HMMRATAC=/mnt/data/userdata/xiangyu/programme/HMMRATAC-1.2.5/HMMRATAC_V1.2.5_exe.jar
java -Xmx30g -jar $HMMRATAC \
-b Merge_ATAC_Arid1a_WT.sorted.bam \
-i Merge_ATAC_Arid1a_WT.sorted.bam.bai \
-g Merge_ATAC_Arid1a_WT.genome.info -o renew_merge_Arid1a_WT_HMMRATAC

HMMRATAC=/mnt/data/userdata/xiangyu/programme/HMMRATAC-1.2.5/HMMRATAC_V1.2.5_exe.jar
java -Xmx30g -jar $HMMRATAC \
-b Merge_ATAC_Arid1a_KO.sorted.bam \
-i Merge_ATAC_Arid1a_KO.sorted.bam.bai \
-g Merge_ATAC_Arid1a_KO.genome.info -o renew_merge_Arid1a_KO_HMMRATAC


library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

files <- list.files(path = "/mnt/data/userdata/abao/project/5_ATAC/1_LHY/NGmerge_produce_fq/filter_bam/filter_MT_bam/final_bam_bw_bai", pattern = "_peaks.gappedPeak$", full.names = TRUE)
name <- gsub("_HMMRATAC_peaks.gappedPeak","",basename(files))

for (i in c(1:length(files))){
  print(i)
  peakAnno <- annotatePeak(files[[i]], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db") 
  peakAnno <- as.data.frame(peakAnno)
  aa <- name[i]
  hh <- paste(aa,"annotatePeak.csv",sep="_")
  bb <- "/mnt/data/userdata/abao/project/5_ATAC/1_LHY/NGmerge_produce_fq/filter_bam/filter_MT_bam/final_bam_bw_bai"
  cc <- paste(bb,hh,sep="/")
  write.csv(peakAnno, file=cc)
}

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)
pdf("Merge_plotAnnoBar.pdf")
plotAnnoBar(peakAnnoList)
dev.off()

pdf("plotDistToTSS.pdf")
plotDistToTSS(peakAnnoList)
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- future_lapply(files, getTagMatrix, windows=promoter)
names(tagMatrixList) <- name

pdf("plotHeatmap.pdf")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

pdf("plotAvgProf.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),free_y=FALSE)
dev.off()


