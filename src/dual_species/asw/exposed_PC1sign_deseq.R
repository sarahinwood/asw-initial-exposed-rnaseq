library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")

asw_dds$group <- factor(paste(asw_dds$Treatment, asw_dds$PC1_sign, sep="_"))
design(asw_dds) <- ~group
asw_dds <- DESeq(asw_dds)

##############
## negative ##
##############

neg_res_group <- results(asw_dds, contrast = c("group", "Exposed_negative", "NC_negative"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
neg_ordered_res_group <- neg_res_group[order(neg_res_group$padj),]
##Make data table and write to output
neg_ordered_res_group_table <- data.table(data.frame(neg_ordered_res_group), keep.rownames = TRUE)
neg_ordered_sig_res_group_table <- subset(neg_ordered_res_group_table, padj < 0.1)

##############
## positive ##
##############

pos_res_group <- results(asw_dds, contrast = c("group", "Exposed_positive", "NC_positive"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
pos_ordered_res_group <- pos_res_group[order(pos_res_group$padj),]
##Make data table and write to output
pos_ordered_res_group_table <- data.table(data.frame(pos_ordered_res_group), keep.rownames = TRUE)
pos_ordered_sig_res_group_table <- subset(pos_ordered_res_group_table, padj < 0.1)

##find 3 DEGs with both comparisons, but different DEGs depending on PC1 sign


plotCounts(asw_dds_location, "", intgroup = c("location"), main="")