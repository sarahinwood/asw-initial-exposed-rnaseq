library("data.table")

sample_data <- fread("data/sample_table.csv", header=TRUE)

##Plotting PCR Target Counts

##ASW##
asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))

##MH##

mh_dds <- readRDS("output/deseq2/mh/mh_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene
plotCounts(mh_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))