library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

sample_data <- fread("data/sample_table.csv", header=TRUE)

#########
## ASW ##
#########

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Treatment))

#PCA
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(asw_vst, intgroup=c("Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=Treatment))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Treatment")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()


##plot counts for PCR target gene
##2 targets in current transcriptome
##ASW_TRINITY_DN4262_c32_g1 - full length
##ASW_TRINITY_DN142_c1_g1 - half as long
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))

########
## MH ##
########

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
mh_dds$group <- factor(paste(mh_dds$abdomen_parasitism_status))

##PCA
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(mh_vst, intgroup=c("abdomen_parasitism_status"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=abdomen_parasitism_status))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Parasitism status")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

##plot counts for PCR target gene
pcr_counts <- plotCounts(mh_dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"), returnData = TRUE)
ggplot(pcr_counts, aes(x=group, y=count))+
  geom_point(size=3, alpha=0.7, color="#440154FF")+
  xlab("Parasitism status")+
  ylab("Normalised count")+
  theme_bw()

##kilA - MH_TRINITY_DN11733_c0_g1
