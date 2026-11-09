# ============================================================================
# Directional N-of-1 MixEnrich -> lnOR / var / z-score matrices
#
# Input : data/cocci_in.rds  -- list(GO_anno, log2fcs)
#           GO_anno : data.frame (path_id, symbol), 5,179 GO BP terms
#           log2fcs : matrix, 7,122 genes x 51 subjects (gene symbols repaired
#                     from Excel date-mangling; e.g. MARCH1, SEPT2, SEP15)
# Output: output/nof1_lnOR_matrix_7122.csv
#         output/nof1_var_matrix_7122.csv
#         output/nof1_zscore_matrix_7122.csv
#           rows = GO:XXXXXXX_up / _down, cols = subjects
#
# Usage : Rscript run_all.R
# ============================================================================
library(mixtools)
library(doParallel)
library(data.table)

source('MixEnrich.R')
source('runMixEnrich.R')
source('differential_expression_direction.R')
source('ef_method_pipeline.R')          # for calc_lnOR_gart

## Load inputs
inputs  <- readRDS('data/cocci_in.rds')
GO_anno <- inputs$GO_anno
log2fcs <- inputs$log2fcs

## N-of-1 MixEnrich: per-subject DEG calls and directional contingency tables
contingency_table_info <- obtain_mixenrich_counts(log2fcs = log2fcs,
                                                  subj_names = colnames(log2fcs),
                                                  poster_prob_cutoff = .99,
                                                  the_alternative = 'two.sided',
                                                  anno_file = GO_anno,
                                                  abs_log2fc_cutoff = log2(1.2),
                                                  num_core = 8)
up_subj_cont_tables   <- contingency_table_info$mixenrich_counts_up
down_subj_cont_tables <- contingency_table_info$mixenrich_counts_down

## Contingency tables -> lnOR / var matrices (rows = GO terms, cols = subjects)
build_ln_mats <- function(count_list) {
  per_subj_lnors <- lapply(count_list, calc_lnOR_gart)
  go_terms <- rownames(per_subj_lnors[[1]])
  per_subj_lnors <- lapply(per_subj_lnors, function(res) res[go_terms, , drop = FALSE])
  lnor_mat <- do.call(cbind, lapply(per_subj_lnors, function(res) res[, 'lnOR', drop = FALSE]))
  var_mat  <- do.call(cbind, lapply(per_subj_lnors, function(res) res[, 'var_lnOR', drop = FALSE]))
  colnames(lnor_mat) <- colnames(var_mat) <- names(count_list)
  rownames(lnor_mat) <- rownames(var_mat) <- go_terms
  list(lnor = lnor_mat, var = var_mat)
}
up_ln   <- build_ln_mats(up_subj_cont_tables)
down_ln <- build_ln_mats(down_subj_cont_tables)
rownames(up_ln$lnor)   <- paste0(rownames(up_ln$lnor), '_up')
rownames(down_ln$lnor) <- paste0(rownames(down_ln$lnor), '_down')
rownames(up_ln$var)    <- paste0(rownames(up_ln$var), '_up')
rownames(down_ln$var)  <- paste0(rownames(down_ln$var), '_down')

lnor_mat   <- rbind(up_ln$lnor, down_ln$lnor)
var_mat    <- rbind(up_ln$var,  down_ln$var)
zscore_mat <- lnor_mat / sqrt(var_mat)

## Write outputs
fwrite(as.data.table(lnor_mat,   keep.rownames = 'GO'), 'output/nof1_lnOR_matrix_7122.csv')
fwrite(as.data.table(var_mat,    keep.rownames = 'GO'), 'output/nof1_var_matrix_7122.csv')
fwrite(as.data.table(zscore_mat, keep.rownames = 'GO'), 'output/nof1_zscore_matrix_7122.csv')
cat('Wrote output/nof1_{lnOR,var,zscore}_matrix_7122.csv\n')
