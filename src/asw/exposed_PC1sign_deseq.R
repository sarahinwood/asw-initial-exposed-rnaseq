library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")

asw_dds$group <- factor(paste(asw_dds$Treatment, asw_dds$PC1_sign, sep="_"))
design(asw_dds) <- ~group
asw_dds <- DESeq(asw_dds)

##############
## negative ##
##############

neg_res_group <- results(asw_dds, contrast = c("group", "Exposed_negative", "NC_negative"), lfcThreshold = 1, alpha = 0.05)
##Order based of padj
neg_ordered_res_group <- neg_res_group[order(neg_res_group$padj),]
##Make data table and write to output
neg_ordered_res_group_table <- data.table(data.frame(neg_ordered_res_group), keep.rownames = TRUE)
neg_ordered_sig_res_group_table <- subset(neg_ordered_res_group_table, padj < 0.05)
neg_sig_annots <- merge(neg_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(neg_sig_annots, "output/deseq2/asw/exposed_pc1_pairwise/PC1_negative_sig_w_annots.csv")

##############
## positive ##
##############

pos_res_group <- results(asw_dds, contrast = c("group", "Exposed_positive", "NC_positive"), lfcThreshold = 1, alpha = 0.05)
##Order based of padj
pos_ordered_res_group <- pos_res_group[order(pos_res_group$padj),]
##Make data table and write to output
pos_ordered_res_group_table <- data.table(data.frame(pos_ordered_res_group), keep.rownames = TRUE)
pos_ordered_sig_res_group_table <- subset(pos_ordered_res_group_table, padj < 0.05)
pos_sig_annots <- merge(pos_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(pos_sig_annots, "output/deseq2/asw/exposed_pc1_pairwise/PC1_positive_sig_w_annots.csv")

saveRDS(asw_dds, "output/deseq2/asw/exposed_pc1_pairwise/asw_dds.rds")

##overlap of DEGs from both analyses
vd <- venn.diagram(x = list("Negative"=neg_sig_annots$rn, "Positive"=pos_sig_annots$rn), filename=NULL, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)

plotCounts(asw_dds_location, "", intgroup = c("location"), main="")