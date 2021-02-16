library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/best_annot_per_gene.csv")

asw_dds$group <- factor(paste(asw_dds$Weevil_Location, asw_dds$PC1_sign, sep="_"))
asw_dds_location <- copy(asw_dds)
design(asw_dds_location) <- ~group
asw_dds_location <- DESeq(asw_dds_location)

##############
## negative ##
##############

neg_res_group <- results(asw_dds_location, contrast = c("group", "Ruakura_negative", "Invermay_negative"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
neg_ordered_res_group <- neg_res_group[order(neg_res_group$padj),]
##Make data table and write to output
neg_ordered_res_group_table <- data.table(data.frame(neg_ordered_res_group), keep.rownames = TRUE)
neg_ordered_sig_res_group_table <- subset(neg_ordered_res_group_table, padj < 0.1)
neg_sig_annots <- merge(neg_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(neg_sig_annots, "output/deseq2/asw/location_pc1_pairwise/PC1_negative_sig_w_annots.csv")

##############
## positive ##
##############

pos_res_group <- results(asw_dds_location, contrast = c("group", "Ruakura_positive", "Invermay_positive"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
pos_ordered_res_group <- pos_res_group[order(pos_res_group$padj),]
##Make data table and write to output
pos_ordered_res_group_table <- data.table(data.frame(pos_ordered_res_group), keep.rownames = TRUE)
pos_ordered_sig_res_group_table <- subset(pos_ordered_res_group_table, padj < 0.1)
pos_sig_annots <- merge(pos_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(pos_sig_annots, "output/deseq2/asw/location_pc1_pairwise/PC1_positive_sig_w_annots.csv")


plotCounts(asw_dds_location, "", intgroup = c("location"), main="")

