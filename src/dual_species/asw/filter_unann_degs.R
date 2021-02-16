#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)

loc_exposed_int_degs <- snakemake@input[['loc_ex_int_degs']]
loc_pairwise_degs <- snakemake@input[['loc_degs']]

loc_exposed_int <- fread("output/deseq2/asw_dual/location_exposure_int/sig_w_annots.csv", na.string="")
loc_pairwise <- fread("output/deseq2/asw_dual/location_pairwise/sig_w_annots.csv", na.string="")
##join two tables together
all_degs <- rbind(loc_pairwise, loc_exposed_int)
all_degs_unann <- all_degs[is.na(sprot_Top_BLASTX_hit),]
fwrite(list(all_degs_unann$rn), snakemake@output[['unann_degs_list']])

# log
sessionInfo()