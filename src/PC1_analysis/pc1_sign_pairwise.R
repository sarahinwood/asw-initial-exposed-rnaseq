library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##filter out genes with less than 10 reads mapping to them - helpful when transcriptome is large
keep <- rowSums(counts(asw_dds)) >= 10
asw_dds_fil <- asw_dds[keep,]

asw_dds_fil$pc1_sign <- factor(paste(asw_dds_fil$PC1_sign))
asw_dds_pc1 <- copy(asw_dds_fil)
design(asw_dds_pc1) <- ~pc1_sign
asw_dds_pc1 <- DESeq(asw_dds_pc1)

res_group <- results(asw_dds_pc1, contrast = c("pc1_sign", "positive", "negative"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw_dual_pc1/res_group.csv")

ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.01)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual_pc1/sig_degs.csv")

dedup_trinotate_report <- fread("/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-transcriptome/output/trinotate/sorted/best_annot_per_gene.csv", na.strings = ".")
ordered_sig_res_group_table$`#gene_id` <- tstrsplit(ordered_sig_res_group_table$rn, "ASW_", keep=c(2))
sig_annots <- merge(ordered_sig_res_group_table, dedup_trinotate_report, by="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual_pc1/sig_annots.csv")

##2107 DEGs
