##
source_path <- '/Users/josoga2/Documents/wale_docs/phd/HB/lncrna_hcc/LNCRNA'
setwd(source_path)

#source base script
source('~/Documents/wale_docs/phd/HB/lncrna_hcc/LNCRNA/deseq2_subset_pipeline.R')

#get lnc annotation
lnc_annotation <- read.delim('mart_export_all_lnc.txt', header = T)

#stage B vs stage c before treatment 
res_stageB_stageC_pre <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "PRE"
  ),
  contrast_col = "BCLC_Stage",
  contrast_levels = c("Stage B", "Stage C"),
  outdir = "cleanup_pass"
)


#stage B vs stage c after treatment 
res_stageB_stageC_post <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST"
  ),
  contrast_col = "BCLC_Stage",
  contrast_levels = c("Stage B", "Stage C"),
  outdir = "cleanup_pass"
)

#responder vs non responder before treatment 
res_resp_nonresp_pre <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "PRE"
  ),
  contrast_col = "Group",
  contrast_levels = c("Responder", "Non-responder"),
  outdir = "cleanup_pass"
)


#responder vs non responder after treatment 
res_resp_nonresp_post <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST"
  ),
  contrast_col = "Group",
  contrast_levels = c("Responder", "Non-responder"),
  outdir = "cleanup_pass"
)

## ==> correlation for this subset
res_resp_nonresp_post_lnccorr <- makeCorrelation(direction = res_resp_nonresp_post$UP, 
                                                 for_corr = res_resp_nonresp_post$FOR_CORR, 
                                                 cutoff = 0.9, 
                                                 n_samp_min = 7)

for (gene in unique(c(res_resp_nonresp_post_lnccorr$gene1, res_resp_nonresp_post_lnccorr$gene2))) {
  cat(gene, sep = '\n ')
}

#==>Chromatin binding, Nucleosome binding



#responder vs non responder after treatment stage B
res_resp_nonresp_post_C <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST",
    BCLC_Stage = "Stage B"
  ),
  contrast_col = "Group",
  contrast_levels = c("Responder", "Non-responder"),
  outdir = "cleanup_pass"
) #==> not so much



#responder vs non responder after treatment stage C
res_resp_nonresp_post_C <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST",
    BCLC_Stage = "Stage C"
  ),
  contrast_col = "Group",
  contrast_levels = c("Responder", "Non-responder"),
  outdir = "cleanup_pass"
)

## ==> correlation for this subset
res_resp_nonresp_post_C_lnccorr_up <- makeCorrelation(direction = res_resp_nonresp_post_C$UP, 
                                                 for_corr = res_resp_nonresp_post_C$FOR_CORR, 
                                                 cutoff = 0.9, 
                                                 n_samp_min = 7)

for (gene in unique(c(res_resp_nonresp_post_C_lnccorr_up$gene1, res_resp_nonresp_post_C_lnccorr_up$gene2))) {
  cat(gene, sep = '\n ')
} #==> mrna binding, protein binding, nucleosome binding, spliceosome, transcription factor activity


res_resp_nonresp_post_C_lnccorr_down <- makeCorrelation(direction = res_resp_nonresp_post_C$DOWN, 
                                                      for_corr = res_resp_nonresp_post_C$FOR_CORR, 
                                                      cutoff = 0.9, 
                                                      n_samp_min = 7)

for (gene in unique(c(res_resp_nonresp_post_C_lnccorr_down$gene1, res_resp_nonresp_post_C_lnccorr_down$gene2))) {
  cat(gene, sep = '\n ')
} #==> mrna binding, protein binding, nucleosome binding, spliceosome, transcription factor activity


#pre vs post in stage b
res_pre_post_responder_B <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Group = 'Responder'
  ),
  contrast_col = "Condition",
  contrast_levels = c("POST", "PRE"),
  outdir = "cleanup_pass"
) #==> correlation for this subset


res_pre_post_responder_B_lnccorr_up <- makeCorrelation(direction = res_pre_post_responder_B$UP, 
                                                        for_corr = res_pre_post_responder_B$FOR_CORR, 
                                                        cutoff = 0.9, 
                                                        n_samp_min = 10)

for (gene in unique(c(res_pre_post_responder_B_lnccorr_up$gene1, res_pre_post_responder_B_lnccorr_up$gene2))) {
  cat(gene, sep = '\n ')
} #==>


