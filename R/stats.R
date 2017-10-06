# 
# ## function to fit the model including both library and guide as parameters
# ## include error handling to try using a different optimizer if the first throws a warning
# #' @import optimx
# #' @importFrom lme4 lmer
# .fitLibGuide <- function(dat) {
#   tryCatch( lmer(log(NormCounts + 1) ~ Treatment + (1 | Library / GuideID), data = as.data.frame(dat)),
#             warning = function(w) { ## if we get a warning, run with a different but slower optimiser
#               suppressMessages(require(optimx))
#               lmer(log(NormCounts + 1) ~ Treatment + (1 | Library / GuideID), data = as.data.frame(dat), 
#                            control = lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
#             }
#   )
# }
# 
# 
# #' fitting the mixed effects model
# #' 
# #' model fitting
# #' 
# #' @return data.table of genes with log fold-change and adjusted p.values
# #' @export
# #' @importFrom broom tidy
# fitModel <- function(countsTable) {
#   
#   modelResults <-
#     filter(countsTable, grep('rg', GuideID, invert = TRUE)) %>%
#     #filter(grep("A1BG", Gene)) %>%
#     #select(GuideID, Treatment, Library, Gene, NormCounts) %>%
#     group_by(Gene) %>% do(
#       if (length(unique(Library)) == 1) {
#         if (length(unique(GuideID) == 1)) {
#           ## if there is only one guide for the gene, return NULL
#           NULL
#         } else {
#           tidy(lmer(log(NormCounts + 1) ~ Treatment + (1 | GuideID), data = .))
#         }
#       } else {
#         tidy( .fitLibGuide(dat = .) )
#       }) %>%
#     ungroup()
#   
#   ## select the rows we want to rename columns for b and t statisics.
#   modelResults <- filter(modelResults, term == "Treatment") %>% 
#     select(Gene, estimate, statistic) %>% 
#     rename(b_stat = estimate, t_stat = statistic)
# 
#   ## here we transform the results to adjusted p-values
#   modelResults <- mutate(modelResults, t_trans = t.transform(t_stat)) %>% 
#     mutate(p_val = p.transform(t_trans)) %>% 
#     mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
#     mutate(fc = exp(b_stat)) %>%
#     arrange(fdr)
#   
#   return(modelResults)
# }
# 
# 
# #' fitting the mixed effects model
# #' 
# #' model fitting
# #' 
# #' @return data.table of genes with log fold-change and adjusted p.values
# #' @export
# #' @importFrom broom tidy
# #' @importFrom lme4 lmer
# fitModel.mc <- function(countsTable, ncores = NA) {

#     require(multidplyr)
#     cluster <- create_cluster(cores = ncores)
#   
#     ct <- ungroup(countsTable) %>% 
#       filter(grep('rg', GuideID, invert = TRUE)) %>%
#         select(GuideID, Treatment, Library, Gene, NormCounts)
#     ct <- data.frame(ct)
# 
#     modelResults <- 
#       partition(ct, Gene, cluster = cluster) %>%
#       do( 
#         if (nrow(unique(.[,'Library'])) == 1) {
#           if (nrow(unique(.[,'GuideID'])) == 1) { 
#             data.frame()
#           } else {
#             broom::tidy(lme4::lmer(log(NormCounts + 1) ~ Treatment + (1 | GuideID), data = .)) 
#           } 
#         } else {
#           broom::tidy( geckoR:::.fitLibGuide(dat = .) )
#         }
#       ) %>%
#       collect() %>%
#       ungroup()
# 
#   ## select the rows we want to rename columns for b and t statisics.
#   modelResults <- filter(modelResults, term == "Treatment") %>% 
#     select(Gene, estimate, statistic) %>% 
#     rename(b_stat = estimate, t_stat = statistic)
#   
#   ## here we transform the results to adjusted p-values
#   modelResults <- mutate(modelResults, t_trans = t.transform(t_stat)) %>% 
#     mutate(p_val = p.transform(t_trans)) %>% 
#     mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
#     mutate(fc = exp(b_stat)) %>%
#     arrange(fdr)
#   
#   return(modelResults)
# }
# 
# #' @importFrom car powerTransform yjPower
# t.transform <- function(x) {
#   lambda <- coef(powerTransform(x, family = "yjPower", method = "BFGS"))
#   tTrans <- yjPower(x, lambda = lambda)
#   return(tTrans)
# }
# 
# #' @importFrom stats pnorm
# p.transform <- function(x) {
#   return(2-2*pnorm(abs(x)))
# }
