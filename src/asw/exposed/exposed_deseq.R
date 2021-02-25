library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")

asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds_exposed <- copy(asw_dds)
design(asw_dds_exposed) <- ~treatment
asw_dds_exposed <- DESeq(asw_dds_exposed)

res_group <- results(asw_dds_exposed, contrast = c("treatment", "Exposed", "NC"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw/exposed/res_group.csv")

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/exposed/sig_annots.csv")
saveRDS(asw_dds_exposed, "output/deseq2/asw/exposed/asw_dds_exposed.rds")

