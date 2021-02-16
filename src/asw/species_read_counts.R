library(data.table)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
counts_table_asw <- (data.table(counts(asw_dds)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("sample_name", "reads_mapped_ASW"))
##merge with csv with total reads
total_reads <- fread("data/total_reads.csv")
total_reads_and_counts <- merge(counts_colSums_asw, total_reads,by="sample_name")
total_reads_and_counts$`mapping_%` <- ((total_reads_and_counts$reads_mapped_ASW/total_reads_and_counts$half_reads_out_bbduk)*100)
fwrite(total_reads_and_counts, "output/deseq2/asw/total_reads_and_counts.csv")

