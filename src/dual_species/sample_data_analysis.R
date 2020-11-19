library("data.table")
library("deseq2")

sample_data <- fread("data/sample_table.csv", header=TRUE)

#########
## ASW ##
#########

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Treatment))

##plot counts for PCR target gene
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))

#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("Cleaned"))

########
## MH ##
########

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Treatment))

##plot counts for PCR target gene
plotCounts(mh_dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))

##PCA
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(mh_vst, intgroup=c("Weevil_Location"))
