library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")

asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds_exposed <- copy(asw_dds)
design(asw_dds_exposed) <- ~treatment
asw_dds_exposed <- DESeq(asw_dds_exposed)

res_group <- results(asw_dds_exposed, contrast = c("treatment", "Exposed", "NC"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)
##nothing siugnificantly differentially expressed

