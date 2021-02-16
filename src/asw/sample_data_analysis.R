library("data.table")
library("DESeq2")

sample_data <- fread("data/sample_table.csv", header=TRUE)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Treatment))

##plot counts for PCR target gene
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))

#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("Weevil_Location"))

