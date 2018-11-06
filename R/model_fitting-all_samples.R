## function to fit the model including both library and guide as parameters
## include error handling to try using a different optimizer if the first throws a warning
#' @import optimx
#' @importFrom lmerTest lmer
#' @importClassesFrom lmerTest lmerModLmerTest
fitPanCellLibGuide <- function(dat) {
  if(sum(dat[,'norm_counts']) < 5) {
    return(data.frame())
  }
  m1 <- tryCatch(lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | library / guide_id), data = as.data.frame(dat), REML = FALSE),
                 warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                   suppressMessages(require(optimx))
                   lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | library / guide_id), data = as.data.frame(dat), REML = FALSE,
                        control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                 }
  )
  
  m1_null <- tryCatch(lmer(log(norm_counts + 1) ~ sample_id + (1 | library / guide_id), data = as.data.frame(dat)),
                      warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                        suppressMessages(require(optimx))
                        lmer(log(norm_counts + 1) ~ sample_id + (1 | library / guide_id), data = as.data.frame(dat), 
                             control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                      }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null)$"Pr(>Chisq)"[2])
  res <- filter(res, term == "treatment")
  res <- mutate(res, anova_pval = anova)
  return(res)
}


#' @import optimx
#' @importFrom lmerTest lmer 
fitPanCellGuide <- function(dat) {
  if(sum(dat[,'norm_counts']) < 5) {
    return(data.frame())
  }
  m1 <- tryCatch(lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | guide_id), data = as.data.frame(dat), REML = FALSE),
                 warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                   suppressMessages(require(optimx))
                   lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | guide_id), data = as.data.frame(dat), REML = FALSE, 
                        control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                 }
  )
  m1_null <- tryCatch(lmer(log(norm_counts + 1) ~  sample_id + (1 | guide_id), data = as.data.frame(dat), REML = FALSE),
                      warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                        suppressMessages(require(optimx))
                        lmer(log(norm_counts + 1) ~  sample_id + (1 | guide_id), data = as.data.frame(dat), REML = FALSE, 
                             control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                      }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null)$"Pr(>Chisq)"[2])
  res <- filter(res, term == "treatment")
  res <- mutate(res, anova_pval = anova)
  return(res)
}

## Before running the model, we have to filter only genes that are tested in at least two cell lines.
# It doesn't make sense to run the model at all for genes that are only tested in one experiment.
# this is especially relevant if different libraries have been used in different experiments!

filterGenesByCellLines <- function(countstable) {
  test <-  group_by(countstable,gene_id) %>%
    summarise(cell_line_counts = n_distinct(cell_line))
  res <-  filter(countstable, gene_id %in% test$gene_id[test$cell_line_counts>1])
  return (res)
}

#' fitting the mixed effects model
#' 
#' model fitting
#' 
#' @return data.table of genes with log fold-change and adjusted p.values
#' @export
PanCell.fitModel <- function(countsTable, controlString = 'control') {
  
  ## As a first step, we have to filter only genes that are tested in at least two cell lines. It doesn't make sense to run the model at all for genes that are only tested in one experiment.
  # This should not influence in cases where all experiments are done with the same library.
  # But will cause problems if we try to integrate Shalem 2014 (GeCKO) with Hart 2015 (TKO)!!
  
  # write code here!!!
  
  
  if(!'norm_counts' %in% colnames(countsTable)) {
    countsTable <- mutate(countsTable, norm_counts = counts)
  }
  
  ct <- ungroup(countsTable) %>% 
    filter(!grepl(controlString, gene_id, ignore.case = TRUE)) %>%
    dplyr::select(guide_id, treatment, library,cell_line, gene_id, norm_counts)
  ct <- data.frame(ct)
  
  modelResults <- ct %>%
    group_by(gene_id) %>%  do(
      if (nrow(unique(.[,'library'])) == 1) {
        if (nrow(unique(.[,'guide_id'])) == 1) {
          ## if there is only one guide for the gene, return NULL
          data.frame()
        } else {
          fitPanCellGuide(dat = .)
        }
      } else {
        ## if we have 2 libraries, but only have one guide per library stick 
        ## with the simpler model
        if( nrow(unique(.[,'guide_id'])) == nrow(unique(.[,'library'])) ) {
          fitPanCellGuide(dat = .)
        } else {
          fitPanCellLibGuide(dat = .)
        }
      }) %>%
    ungroup()
  
  # select the rows we want to rename columns for b and t statisics.
  modelResults <- select(modelResults, gene_id, Estimate, t.value, anova_pval) %>%
    rename(b_stat = Estimate, t_stat = t.value) %>%
    mutate(fdr_anova = p.adjust(anova_pval, method = "fdr")) %>%
    mutate(fc = exp(b_stat)) %>%
    arrange(fdr_anova)
  
  return(modelResults)
}