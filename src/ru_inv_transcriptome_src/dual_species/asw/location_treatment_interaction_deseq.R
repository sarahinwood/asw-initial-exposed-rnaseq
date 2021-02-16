library(tximport)
library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")

##design
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$treatment <- factor(paste(asw_dds$Treatment))
asw_dds$clean <- factor(paste(asw_dds$Cleaned))
design(asw_dds) <- ~clean+location+treatment+location:treatment

##want interaction between location:treatment
##WT as all factors only have 2 levels - WT is for 3+
dds_WT_asw <- DESeq(asw_dds, test="Wald")
resultsNames(dds_WT_asw)
res_group <- results(dds_WT_asw, alpha=0.1)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)

trinotate <- fread("data/asw-transcriptome/output/trinotate/sorted/best_annot_per_gene.csv")
trinotate$edited_id <- paste("ASW", trinotate$`#gene_id`, sep="_")
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="edited_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/location_exposure_int/sig_w_annots.csv")

plotCounts(asw_dds_location, "", intgroup = c("location"), main="")

##plot counts for genes of interest, sub in name
plotCounts(dds_WT_asw, "ASW_TRINITY_DN3182_c0_g1", intgroup = c("location", "treatment"))
