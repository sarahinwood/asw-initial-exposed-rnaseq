library(tximport)
library(data.table)
library(DESeq2)

mh_gene_trans_map <- "data/asw-mh-combined-transcriptome/output/mh_edited_transcript_ids/Trinity.fasta.gene_trans_map"

gene2tx <- fread(mh_gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
mh_dds <- "output/deseq2/mh_dual/mh_dual_dds.rds"
saveRDS(dds, mh_dds)
