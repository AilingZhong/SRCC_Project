project <- c("Tumor")
deal_design <- c("TPA","TP")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- paste(getwd(),"/0_DESeq2_Tumor_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/ebg_mm10.RData")
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/txdb_mm10.RData")
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
csvfile <- file.path(indir, "Tumor_sampletable.csv")
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
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","TPA","TP"))
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
names(all_summry) <- c("TPA_1","TPA_2","TPA_3","TP_1","TP_2","DESeq2_TPA_1","DESeq2_TPA_2","DESeq2_TPA_3","DESeq2_TP_1","DESeq2_TP_2", "baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","entrez","GENENAME")
write.csv(all_summry, "renew_DEseq2normalized_TPA_VS_TP_allsummry.csv")



all_summry <- cbind(countdata,res_1)
all_summry <- na.omit(all_summry)
p005 <- subset(all_summry,pvalue < 0.05)
p005 <- subset(p005,log2FoldChange < -1 | log2FoldChange > 1 )
p005 <- p005[order(p005[,14]),]
names(p005) <- c("TPA_1","TPA_2","TPA_3","TP_1","TP_2","DESeq2_TPA_1","DESeq2_TPA_2","DESeq2_TPA_3","DESeq2_TP_1","DESeq2_TP_2", "baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","entrez","GENENAME")
write.csv(p005,"renew_TPA_VS_TP_p005_log1_forheatmap.csv")



# first_1
res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),9]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),6]
upres_1 <- res_1[with(res_1,y>=0.8),]
downres_1  <- res_1[with(res_1,y<= -0.8),]


ee	<-as.matrix(upres_1$entrez)
dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, organism = organism, keyType = "ncbi-geneid",pvalueCutoff = 0.05,pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,qvalueCutoff = 0.5, use_internal_data = FALSE)
KEGGupres_1 <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")

ee	<-as.matrix(downres_1$entrez)
dd <- as.vector(ee)
KEGGdownres_1 <- enrichKEGG(gene =dd, organism = organism, keyType = "ncbi-geneid",pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = 10, maxGSSize = 500,qvalueCutoff = 0.5, use_internal_data = FALSE)
KEGGdownres_1 <- setReadable(KEGGdownres_1, "org.Mm.eg.db", keyType="ENTREZID")


write.csv(KEGGupres_1, file = KEGGupres_1_file.csv)
write.csv(KEGGdownres_1, file = KEGGdownres_1_file.csv)
pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file =GOres_1_all_UP_pdf)
ff1
dev.off()

ee <- as.matrix(downres_1$entrez)
dd <- as.vector(ee)
GOdownres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)
ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file =GOres_1_all_DOWN_pdf)
ff1
dev.off()





TPA_VS_TP_count_tpm <- read.csv("/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/Tumor_workfile/Tumor_5_TPA_VS_TP_count_tpm_symbol_and_anno.csv")

convert_mouse_to_human_GeneList <- function(x){
	require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), 
    martL = mouse, uniqueRows=T)
  print(head(genesV2))
  return(genesV2)}
mouse_to_human <- convert_mouse_to_human_GeneList(TPA_VS_TP_count_tpm$X)

mousetohuman <- merge(mouse_to_human,TPA_VS_TP_count_tpm,by.x="MGI.symbol",by.y="X")
write.csv(mousetohuman,"5_tumor_orgenoid_mousetohuman.csv")


GO_DOWN_geneID <- read.csv("/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/Tumor_workfile/GO_DOWN_geneID.csv")
TPA_VS_TP_count_tpm_symbol_and_anno <- read.csv("/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/Tumor_workfile/Tumor_group_KEGGandGO/Tumor_5_TPA_VS_TP_count_tpm_symbol_and_anno.csv")
aa <- merge(TPA_VS_TP_count_tpm_symbol_and_anno,GO_DOWN_geneID,by.x="X",by.y="GO_DOWN_geneID")
write.csv(aa,"5_tumor_GO_down_genesets.csv")




/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct
/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls



java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_2vs2_Normal_workfile/vesicle.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/tumor_vesicle -gui false


java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_2vs2_Normal_workfile/exocytosis.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/tumor_exocytosis -gui false


java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_2vs2_Normal_workfile/Golgi.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/tumor_Golgi -gui false






/mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/1_Old_version_TPM/Tumor_workfile/Tumor_GSEA/CALCIUM.gmt



java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/1_Old_version_TPM/Tumor_workfile/Tumor_GSEA/CALCIUM.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/tumor_CALCIUM -gui false





java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Mouse_SRCC_EMT.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/tumor_Metastasis -gui false





java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_DEseq2normalized_TPA_VS_TP.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/Tumor_Arid1a.cls#TPA_versus_TP \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/DNA_packaging.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/DNA_packaging -gui false



java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/reference_data/TCGA/SRCC_VS_Intest_GSEA/SRCC12_Aden74_for.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/reference_data/TCGA/SRCC_VS_Intest_GSEA/SRCC12_Aden74_for.cls#SRCC_versus_Intest \
-gmx /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/DNA_packaging.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 1000 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/BAM_MM10/0_DESeq2_Tumor_workfile/DESeq2_Tumor_GSEA/SRCC_versus_Intest_DNA_packaging -gui false




















