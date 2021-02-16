library("data.table")
library("DESeq2")

sample_data <- fread("sample_table_pca_info.csv", header=TRUE)
sample_data$extraction_batch <- tstrsplit(sample_data$sample_name, "_", keep=c(3))

asw_dds <- readRDS("output/inv-ru/deseq2/asw/asw_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Cleaned))

#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("group"))
pca_table <- plotPCA(asw_vst, intgroup=c("group"), returnData=TRUE)
fwrite(pca_table, "output/inv-ru/deseq2/asw/pca_analysis/pca_table.csv")


##plot counts for PCR target gene
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))