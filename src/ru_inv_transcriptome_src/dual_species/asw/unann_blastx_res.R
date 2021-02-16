library("data.table")
library("dplyr")

asw_nr_blastx <- fread("output/deseq2/asw_dual/unann/nr_blastx.outfmt3")
setnames(asw_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(asw_nr_blastx, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
asw_min_evalues <- asw_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
asw_min_evalues$rn <- tstrsplit(asw_min_evalues$transcript_id, "_i", keep=c(1))

##merge with DEG lists
loc_ex_int_degs <- fread('output/deseq2/asw_dual/location_exposure_int/sig_w_annots.csv')
loc_ex_blast <- merge(loc_ex_int_degs, asw_min_evalues, by="rn", all.x=TRUE)
fwrite(loc_ex_blast, 'output/deseq2/asw_dual/location_exposure_int/sig_w_blast_annots.csv')

loc_degs <- fread('output/deseq2/asw_dual/location_pairwise/sig_w_annots.csv')
loc_blast <- merge(loc_degs, asw_min_evalues, by="rn", all.x=TRUE)
fwrite(loc_blast, "output/deseq2/asw_dual/location_pairwise/sig_w_blast_annots.csv")

