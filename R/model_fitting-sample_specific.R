## function to fit the model including both library and guide as parameters
## include error handling to try using a different optimizer if the first throws a warning
#' @import optimx
#' @importFrom lmerTest lmer
#' @importClassesFrom lmerTest lmerModLmerTest
.fitLibGuide.sampleSpecific <- function(dat) {
  if(sum(dat[,'norm_counts']) < 5) {
    return(data.frame())
  }
  m1 <- tryCatch(
    suppressMessages(
      lmer(log(norm_counts + 1) ~ treatment + sample_id + treatment:sample_id + (1 | library / guide_id), 
           data = as.data.frame(dat), REML = FALSE)
    ),
    warning = function(w) { ## if we get a warning, run with a different but slower optimiser
      suppressMessages(require(optimx))
      lmer(log(norm_counts + 1) ~ treatment + sample_id + treatment:sample_id + (1 | library / guide_id), 
           data = as.data.frame(dat), REML = FALSE,
           control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
    }
  )
  
  m1_null <- tryCatch(
    suppressMessages(
      lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | library / guide_id), 
           data = as.data.frame(dat), REML = FALSE)
    ),
    warning = function(w) { ## if we get a warning, run with a different but slower optimiser
      suppressMessages(require(optimx))
      lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | library / guide_id), 
           data = as.data.frame(dat), REML = FALSE,
           control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
    }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null,refit=FALSE)$"Pr(>Chisq)"[2])
  res <- mutate(res, anova_pval = anova)
  return(res)
}


#' @import optimx
#' @importFrom lmerTest lmer 
.fitGuide.sampleSpecific <- function(dat) {
  if(sum(dat[,'norm_counts']) < 5) {
    return(data.frame())
  }
  m1 <- tryCatch(
    suppressMessages(
      lmer(log(norm_counts + 1) ~ treatment + sample_id + treatment:sample_id + (1 | guide_id), 
           data = as.data.frame(dat), REML = FALSE)
    ),
    warning = function(w) { ## if we get a warning, run with a different but slower optimiser
      suppressMessages(require(optimx))
      lmer(log(norm_counts + 1) ~ treatment + sample_id + treatment:sample_id + (1 | guide_id), 
           data = as.data.frame(dat), REML = FALSE, 
           control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
    }
  )
  m1_null <- tryCatch(
    suppressMessages(
      lmer(log(norm_counts + 1) ~ treatment + sample_id + (1 | guide_id), 
           data = as.data.frame(dat), REML = FALSE)
    ),
    warning = function(w) { ## if we get a warning, run with a different but slower optimiser
      suppressMessages(require(optimx))
      lmer(log(norm_counts + 1) ~  treatment + sample_id + (1 | guide_id), 
           data = as.data.frame(dat), REML = FALSE, 
           control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
    }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null, refit=FALSE)$"Pr(>Chisq)"[2])
  res <- mutate(res, anova_pval = anova)
  return(res)
}

#' fitting the mixed effects model
#' 
#' model fitting
#' 
#' @return data.table of genes with log fold-change and adjusted p.values
#' @export
memcrispr.fitModel.sampleSpecific <- function(countsTable, controlString = 'control') {
  
  if(!'norm_counts' %in% colnames(countsTable)) {
    countsTable <- mutate(countsTable, norm_counts = counts)
  }
  
  ct <- ungroup(countsTable) %>% 
    filter(!grepl(controlString, gene_id, ignore.case = TRUE)) %>%
    select(guide_id, treatment, library, sample_id, gene_id, norm_counts)
  ct <- data.frame(ct)
  
  modelResults <- ct %>%
    group_by(gene_id) %>%  do(
      if (nrow(unique(.[,'library'])) == 1) {
        if (nrow(unique(.[,'guide_id'])) == 1) {
          ## if there is only one guide for the gene, return NULL
          data.frame()
        } else {
          MEMcrispR:::.fitGuide(dat = .)
        }
      } else {
        ## if we have 2 libraries, but only have one guide per library stick 
        ## with the simpler model
        if( nrow(unique(.[,'guide_id'])) == nrow(unique(.[,'library'])) ) {
          MEMcrispR:::.fitGuide.sampleSpecific(dat = .)
        } else {
          MEMcrispR:::.fitLibGuide.sampleSpecific(dat = .)
        }
      }) %>%
    ungroup()
  
  ## select the rows we want to rename columns for b and t statisics.
  modelResults <- select(modelResults, gene_id, term, Estimate, t.value, anova_pval) %>% 
    rename(b_stat = Estimate, t_stat = t.value) %>% 
    mutate(fdr_anova = p.adjust(anova_pval, method = "fdr")) %>%
    mutate(fc = exp(b_stat)) %>%
    arrange(gene_id)
  
  return(modelResults)
}
