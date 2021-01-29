
## function to fit the model including both library and guide as parameters
## include error handling to try using a different optimizer if the first throws a warning
#' @import optimx
#' @importFrom lmerTest lmer
#' @importClassesFrom lmerTest lmerModLmerTest
.fitSampleLibGuide <- function(dat) {
  if(sum(dat[,'norm_counts']) == 0) {
    return(data.frame())
  }
  m1 <- tryCatch( summary(lmer(log(norm_counts + 1) ~ treatment + (1 | library / guide_id), 
                                data = as.data.frame(dat))),
                   warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                     suppressMessages(require(optimx))
                     summary(lmer(log(norm_counts + 1) ~ treatment + (1 | library / guide_id), 
                                  data = as.data.frame(dat), 
                                  control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B"))))
                   }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null)$"Pr(>Chisq)"[2])
  res <- mutate(res, anova_pval = anova)
  return(res)
}




## function to fit the model including both library and guide as parameters
## include error handling to try using a different optimizer if the first throws a warning
#' @import optimx
#' @importFrom lmerTest lmer
#' @importClassesFrom lmerTest lmerModLmerTest
.fitLibGuide <- function(dat) {
  if(sum(dat[,'norm_counts']) == 0) {
    return(data.frame())
  }
  m1 <- tryCatch( lmerTest::lmer(log(norm_counts + 1) ~ treatment + (1 | library / guide_id), 
                                data = as.data.frame(dat), REML = FALSE),
            warning = function(w) { ## if we get a warning, run with a different but slower optimiser
              suppressMessages(require(optimx))
              lmerTest::lmer(log(norm_counts + 1) ~ treatment + (1 | library / guide_id), data = as.data.frame(dat), 
                           control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
            }
  )
  m1_null <- tryCatch(lmer(log(norm_counts + 1) ~ (1 | library / guide_id), 
                           data = as.data.frame(dat), REML = FALSE),
                      warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                        suppressMessages(require(optimx))
                        lmer(log(norm_counts + 1) ~ (1 | library / guide_id), 
                             data = as.data.frame(dat), REML = FALSE, 
                             control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                      }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <-as.numeric(anova(m1, m1_null)$"Pr(>Chisq)"[2])
  res <- mutate(res, anova_pval = anova)
  return(res)
}


#' @import optimx
#' @importFrom lmerTest lmer 
.fitGuide <- function(dat) {
  if(sum(dat[,'norm_counts']) < 10) {
    return(data.frame())
  }
  m1 <- tryCatch( lmer(log(norm_counts + 1) ~ treatment + (1 | guide_id), 
                          data = as.data.frame(dat), REML = FALSE),
                   warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                     suppressMessages(require(optimx))
                     lmer(log(norm_counts + 1) ~ treatment + (1 | guide_id), 
                          data = as.data.frame(dat), REML = FALSE,
                          control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                   }
  )
  m1_null <- tryCatch(lmer(log(norm_counts + 1) ~ (1 | guide_id), 
                           data = as.data.frame(dat), REML = FALSE),
                      warning = function(w) { ## if we get a warning, run with a different but slower optimiser
                        suppressMessages(require(optimx))
                        lmer(log(norm_counts + 1) ~ (1 | guide_id), 
                             data = as.data.frame(dat), REML = FALSE, 
                             control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
                      }
  )
  res <- coefficients(summary(m1)) %>%
    data.frame(term = rownames(.))
  anova <- as.numeric(anova(m1, m1_null)$"Pr(>Chisq)"[2])
  res <- mutate(res, anova_pval = anova)
  return(res)
}

#' fitting the mixed effects model
#' 
#' model fitting
#' 
#' @return data.table of genes with log fold-change and adjusted p.values
#' @export
memcrispr.fitModel <- function(countsTable, controlString = 'control') {
  
  if(!'norm_counts' %in% colnames(countsTable)) {
    countsTable <- mutate(countsTable, norm_counts = counts)
  }
  
  ct <- ungroup(countsTable) %>% 
    filter(!grepl(controlString, gene_id, ignore.case = TRUE)) %>%
    select(guide_id, treatment, library, gene_id, norm_counts)
  ct <- data.frame(ct)
  
  modelResults <- ct %>%
    group_by(gene_id) %>%  do(
      if (nrow(unique(.[,'library'])) == 1) {
        if (nrow(unique(.[,'guide_id'])) == 1) {
          ## if there is only one guide for the gene, return NULL
          data.frame()
        } else {
          suppressMessages(
            MEMcrispR:::.fitGuide(dat = .)
          )
        }
      } else {
        ## if we have 2 libraries, but only have one guide per library stick 
        ## with the simpler model
        if( nrow(unique(.[,'guide_id'])) == nrow(unique(.[,'library'])) ) {
          suppressMessages(
            MEMcrispR:::.fitGuide(dat = .)
          )
        } else {
          suppressMessages(
            MEMcrispR:::.fitLibGuide(dat = .)
          )
        }
      }) %>%
    ungroup()
  
  ## select the rows we want to rename columns for b and t statisics.
  modelResults <- filter(modelResults, term == "treatment") %>% 
    select(gene_id, Estimate, t.value, Pr...t..) %>% 
    rename(b_stat = Estimate, t_stat = t.value, p_val = Pr...t..) %>% 
    mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
    mutate(fc = exp(b_stat)) %>%
    arrange(fdr)
  
  return(modelResults)
}


#' fitting the mixed effects model
#' 
#' model fitting
#' 
#' @return data.table of genes with log fold-change and adjusted p.values
#' @export
memcrispr.fitModel.mc <- function(countsTable, ncores) {
  
  if(missing(ncores)) {
    ncores <- getOption(Ncpus, 1L)
  }
  
  require(multidplyr)
  cluster <- multidplyr::new_cluster(n = ncores)
  
  if(!'norm_counts' %in% colnames(countsTable)) {
    countsTable <- mutate(countsTable, norm_counts = counts)
  }
  
  ct <- ungroup(countsTable) %>% 
    filter(!grepl('control', gene_id, ignore.case = TRUE)) %>%
    select(guide_id, treatment, library, gene_id, norm_counts)
  ct <- data.frame(ct) %>%
    group_by(gene_id)
  
  modelResults <- 
    partition(ct, cluster = cluster) %>%
    do( 
      if (nrow(unique(.[,'library'])) == 1) {
        if (nrow(unique(.[,'guide_id'])) == 1) { 
          data.frame()
        } else {
          MEMcrispR:::.fitGuide(dat = .) 
        } 
      } else {
        ## if we have 2 libraries, but only have one guide per library stick 
        ## with the simpler model
        if( nrow(unique(.[,'guide_id'])) == nrow(unique(.[,'library'])) ) {
          MEMcrispR:::.fitGuide(dat = .)
        } else {
          MEMcrispR:::.fitLibGuide(dat = .)
        }
      }
    ) %>%
    collect() %>%
    ungroup()
  
  ## select the rows we want to rename columns for b and t statisics.
  modelResults <- filter(modelResults, term == "treatment") %>% 
    select(gene_id, Estimate, t.value, Pr...t..) %>% 
    rename(b_stat = Estimate, t_stat = t.value, p_val = Pr...t..) %>% 
    mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
    mutate(fc = exp(b_stat)) %>%
    arrange(desc(fc))
  
  return(modelResults)
}


