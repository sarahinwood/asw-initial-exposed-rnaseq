library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/best_annot_per_gene.csv")

##design
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$clean <- factor(paste(asw_dds$Cleaned))
design(asw_dds) <- ~clean+location+treatment+location:treatment

##want interaction between location:treatment
##WT as all factors only have 2 levels - LRT is for 3+
dds_WT_asw <- DESeq(asw_dds, test="Wald")
resultsNames(dds_WT_asw)
res_group <- results(dds_WT_asw, alpha=0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

sig_annots <- merge(trinotate, ordered_sig_res_group_table, by.x="#gene_id", by.y="rn", all.y=TRUE)
fwrite(sig_annots, "output/deseq2/asw/location_exposure_int/sig_w_annots.csv")
saveRDS(dds_WT_asw, "output/deseq2/asw/location_exposure_int/dds.rds")


##plot counts for genes of interest, sub in name
plotCounts(dds_WT_asw, "TRINITY_DN9423_c0_g1", intgroup = c("location", "treatment"))
