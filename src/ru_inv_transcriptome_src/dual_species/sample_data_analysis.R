library("data.table")
library("deseq2")

sample_data <- fread("data/sample_table.csv", header=TRUE)

#########
## ASW ##
#########

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Treatment))

##plot counts for PCR target gene
##2 targets in current transcriptome
##ASW_TRINITY_DN4262_c32_g1 - full length
##ASW_TRINITY_DN142_c1_g1 - half as long
plotCounts(asw_dds, "ASW_TRINITY_DN142_c1_g1", intgroup = c("group"))

#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(asw_vst, intgroup=c("Weevil_Location"))

########
## MH ##
########

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Treatment))

##plot counts for PCR target gene
plotCounts(mh_dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))
##kilA - MH_TRINITY_DN11733_c0_g1

##PCA
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(mh_vst, intgroup=c("Weevil_Location"))
