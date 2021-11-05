library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
counts_table_asw <- (data.table(counts(asw_dds)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("sample_name", "reads_mapped_ASW"))
##merge with csv with total reads
total_reads <- fread("data/total_reads.csv")
total_reads_and_counts <- merge(counts_colSums_asw, total_reads,by="sample_name")

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
counts_table_mh <- (data.table(counts(mh_dds)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("sample_name", "reads_mapped_MH"))
##merge with existing total reads and counts table
total_reads_and_counts <- merge(counts_colSums_mh, total_reads_and_counts,by="sample_name")
##calc. total mapped reads
total_reads_and_counts$total_mapped_reads <- (total_reads_and_counts$reads_mapped_MH + total_reads_and_counts$reads_mapped_ASW)
total_reads_and_counts$`mapping_%` <- ((total_reads_and_counts$total_mapped_reads/total_reads_and_counts$half_reads_out_bbduk)*100)
total_reads_and_counts$`%_ASW` <- ((total_reads_and_counts$reads_mapped_ASW/total_reads_and_counts$total_mapped_reads)*100)
total_reads_and_counts$`%_MH` <- ((total_reads_and_counts$reads_mapped_MH/total_reads_and_counts$total_mapped_reads)*100)
total_reads_and_counts$`mh/total_reads` <- ((total_reads_and_counts$reads_mapped_MH/total_reads_and_counts$half_reads_out_bbduk)*100)
total_reads_and_counts$`asw/total_reads` <- ((total_reads_and_counts$reads_mapped_ASW/total_reads_and_counts$half_reads_out_bbduk)*100)
fwrite(total_reads_and_counts, "output/deseq2/asw_dual/total_reads_and_counts.csv")

ggplot(total_reads_and_counts, aes)

