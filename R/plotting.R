# #' @importFrom tidyr spread
# #' @importFrom gridExtra grid.arrange
# #' @importFrom  utils combn
# memcrispr.compareControlGuides <- function(countTable, controlString = NULL, treatment = NULL) {
#   
#   dat <- ungroup(countTable)
#   if(!is.null(controlString)) {
#     dat <- filter(dat, grepl(controlString, gene_id, ignore.case = TRUE))
#   }
#   
#   if(!'norm_counts' %in% colnames(dat)) {
#     dat <- mutate(dat, norm_counts = log10(counts+1))
#   } else {
#     dat <- mutate(dat, norm_counts = log10(norm_counts+1))
#   }
#   
#   rep_list <- split(dat, 
#                     paste(dat[['treatment']], dat[['replicate']], dat[['library']]))
#   pairings <- combn(names(rep_list), 2)
#   
#   tmp <- list()
#   for(i in 1:ncol(pairings)) {
#     tmp[[ i ]] <- full_join(rep_list[[ as.character(pairings[1,i]) ]], 
#                             select(rep_list[[ as.character(pairings[2,i]) ]], guide_id, treatment, library, replicate, norm_counts), 
#                             by = c('guide_id')) %>%
#       mutate(replicates = paste(pairings[1,i], "vs", pairings[2,i])) %>%
#       mutate(comparison = paste0("Treatment ", treatment.x, ", Rep ", replicate.x, "\nvs\n",
#                                  "Treatment ", treatment.y, ", Rep ", replicate.y))
#   }
#   dat <- do.call("rbind", tmp)
#   
#   p1 <- ggplot(dat, aes(x = norm_counts.x, y = norm_counts.y)) +
#     geom_point(colour = "#296E62", alpha = 0.5) +
#     geom_abline(intercept = 0, slope = 1, col = "#A63E4A") +
#     facet_grid(. ~ comparison)
# 
#   ## force the plots to be square
#   p1 <- p1 + 
#     xlim(0, max(p1$data$norm_counts.x, p1$data$norm_counts.y, na.rm = TRUE)) + 
#     ylim(0, max(p1$data$norm_counts.x, p1$data$norm_counts.y, na.rm = TRUE)) + 
#     coord_fixed(ratio = 1) +
#     xlab('log'[10]~' counts') +
#     ylab('log'[10]~' counts') +
#     theme_bw()
#    
#   return(p1)
# }

#' @importFrom GGally lowertriangle uppertriangle
.scatterMatrix <- function(data) {
  
  ## calculate correlation values for plotting in
  ## upper triangle
  a <- uppertriangle(data, corMethod = "pearson")
  a[,'xvalue'] <- a[,'yvalue'] <- max(data)/2
  
  ltdata.new <- lowertriangle(data)
  
  ## create a data frame so we only draw abline in the lower
  ## triangle of plots
  abline_table <- ltdata.new %>% 
    filter(!is.na(xvalue) & !is.na(yvalue)) %>% 
    select(xlab, ylab) %>% 
    distinct() %>% 
    mutate(slope = 1, intercept = 0)
  
  r <- ggplot(ltdata.new, 
              mapping = aes_string(x = "xvalue", y = "yvalue")) + 
    geom_point(colour = "grey20", alpha = 0.5, na.rm = TRUE) +
    geom_abline(aes(slope = slope, intercept = intercept), data = abline_table) +
    geom_text(data = a, aes_string(label = "r"), colour = "black") +
    facet_grid(ylab ~ xlab, scales = "fixed") + 
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlim(0, max(data, na.rm = TRUE)) +
    ylim(0, max(data, na.rm = TRUE))
  
  densities <- do.call("rbind", 
                       lapply(1:ncol(data), 
                              function(i) {
                                data.frame(xlab = names(data)[i], 
                                           ylab = names(data)[i], 
                                           x = data[, i])
                              }))
  for (m in 1:ncol(data)) {
    j <- subset(densities, xlab == names(data)[m])
    r <- r + stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
                          data = j, position = "identity", 
                          geom = "line", color = "black")
  }
  r
}

#' Compare the counts of control guides for all samples in the experiment
#' 
#' @param countTable Standard \code{memcrispr} table of count data.
#' @param controlString Specify a string that identifys the control guides.  
#' This is often 'control', 'non-targetting' etc.  If this is left as NULL
#' all guides will be compared, which can be very slow and/or crash.
#' 
#' @export
#' @importFrom tidyr spread
#' @importFrom GGally ggscatmat
memcrispr.compareControlGuides <- function(countTable, controlString = NULL) {
  
  dat <- ungroup(countTable)
  if(!is.null(controlString)) {
    dat <- filter(dat, grepl(controlString, gene_id, ignore.case = TRUE))
  }
  
  if(!'norm_counts' %in% colnames(dat)) {
    dat <- mutate(dat, norm_counts = log10(counts+1))
  } else {
    dat <- mutate(dat, norm_counts = log10(norm_counts+1))
  }
  
  lib_list <- split(dat, dat[[ 'library' ]])
  plots <- list()
  
  for(i in seq_along(lib_list)) {
    dat2 <- lib_list[[i]] %>% 
      mutate(comparison = paste0("Treatment: ", treatment, "\nRep: ", replicate)) %>%
      select(guide_id, norm_counts, comparison) %>%
      spread(key = comparison, value = norm_counts)
    
    plots[[i]] <- .scatterMatrix(as.data.frame(dat2)[,-1]) +
      xlab('log'[10]~' counts') +
      ylab('log'[10]~' counts') +
      ggtitle(paste("Library:", names(lib_list)[i]))
  }
  return(plots)
}
  
  
#' Plot distribution of guide counts
#' 
#' 
#' 
#' @export
#' @importFrom tidyr gather
memcrispr.guideDistribution <- function(countTable) {
  
  dat <- 
    ungroup(countTable) %>%
    gather(count_type, Counts, which(colnames(countTable) %in% c("counts", "norm_counts")))
  
  ggplot(data = dat, aes(group = interaction(treatment, library, replicate),
                         colour = as.factor(replicate),
                         linetype = library,
                         x = log10(Counts+1))) + 
    scale_colour_discrete('replicate') +
    scale_linetype_discrete('library') +
    geom_density() +
    facet_wrap(~ count_type + treatment, labeller = "label_both") +
    xlab('log'[10]~' counts') +
    theme_bw()

}

#' Create volcano plot comparing fold-change and p-value for each gene
#' 
#' For every gene in a model results table, this function plots the \code{log2}
#' fold change against the \code{-log10} p-value.  Genes that fall below 
#' specified p-value threshold and above a fold change threshold are 
#' highlighted and named in the plot.
#' 
#' @param topTable Standard \code{memcrispr} table of model results.  Typically
#' produced by \code{\link{memcrispr.fitModel}}.
#' @param fc.thresh Specify a fold change cutoff. Expects to be provide the raw
#' fold change, which will be automatically converted to log2 scale for 
#' plotting.
#' @param p.thresh Specify a p-value cutoff.  This exepcts to be provided with
#' the actual p-value e.g. 0.05, which will automatically be converted to 
#' -log10 scale for plotting.
#' 
#' @export
memcrispr.volcanoPlot <- function(topTable, fc.thresh = 1, p.thresh = 0.05) {
  
  ggplot(topTable) + 
    geom_point(aes(x = log2(fc), y = -log10(1e-16+p_val)), alpha = 0.4, colour = "#939393") +
    geom_point(data = filter(topTable, log2(fc) > fc.thresh, p_val < p.thresh),
            aes(x =  log2(fc), y = -log10(1e-16+p_val)), color = "#EA0303", size = 1) +
    geom_point(data = filter(topTable, log2(fc) > fc.thresh, p_val < p.thresh),
               aes(x =  log2(fc), y = -log10(1e-16+p_val)), color = "#EA0303", size = 4, pch = 1) +
    geom_text(data = filter(topTable, log2(fc) > fc.thresh, p_val < p.thresh),
            aes(x =  log2(fc), y = -log10(1e-16+p_val), label = as.character(gene_id)), 
            color = "#000000", size = 3) +
    ## include dotted lines for the thresholds
    geom_hline(yintercept = -log10(p.thresh), linetype = 2) +
    geom_vline(xintercept = fc.thresh, linetype = 2) +
    ## axis labels
    xlab('log'[2]~' fold-change') +
    ylab('-log'[10]~ ' p-value') +
    xlim(-6,6) + 
    theme_bw()

}

#' @export
memcrispr.libraryCounts <- function(countTable) {
  
  dat <- ungroup(countTable) %>%
    group_by(treatment, library, replicate) %>%
    summarize(totalCount = sum(counts)) %>%
    mutate(group = paste0('replicate: ', replicate, '\nlibrary: ', library)) 
  
  
  ggplot(dat, aes(x = group, y = totalCount, fill = as.factor(library))) + 
    geom_bar(stat = "identity") +
    # facet_wrap(~ treatment, labeller = "label_both") +
    ylab("Total Guide Count") +
    xlab("CRISPR Library") + 
    scale_fill_discrete('library') +
    facet_grid(. ~ treatment, scales = "free_x", 
               space = "free_x", labeller = "label_both") +
    theme_bw()
  
}

## internal function to compute the PCA. Since we use this on
## raw and normalized counts sepately we keep it as a distinct function.
#' @importFrom tidyr unite
.pcaInternal <- function(dat, countType = "counts") {
  
  if(countType == "counts") 
    dat <- mutate(dat, counts = log(counts+1))
  else if(countType == "norm_counts")
    dat <- mutate(dat, counts = log(norm_counts+1))     
  else 
    stop("Unspecified counts column")
  
  dat <- dat %>%
    select(treatment, replicate, library, seq, counts) %>% 
    unite(all_factors, treatment, replicate, library, remove = FALSE) %>%
    select(-treatment, -replicate) %>%
    group_by(all_factors) %>%
    distinct(seq, .keep_all = TRUE) %>%
    spread(key = all_factors, value = counts)
  
  d <- NULL
  
  for(lib in unique(dat[['library']])) {
    pca_input <- dat %>% 
      select(grep(pattern = lib, names(.))) %>%
      na.omit() %>%
      as.data.frame()
    
    pca.res <- prcomp(t(pca_input))
    percentVar <- 100 * pca.res$sdev^2 / sum( pca.res$sdev^2 )
    
    d <- rbind(d, data.frame(PC1=pca.res$x[,1], PC2=pca.res$x[,2], 
                             PC1.var = round(percentVar[1], 2), PC2.var = round(percentVar[2], 2),
                             do.call(rbind, strsplit(rownames(pca.res$x), "_")), count_type = countType))
  }
  
  colnames(d)[5:7] <- c("treatment", "replicate", "library")
  d <- mutate(d, group = paste0('treatment: ', treatment, '\nlibrary: ', library))
  
  return(d)
  
}


#' @export
#' @importFrom tidyr unite
memcrispr.pcaPlot <- function(countTable, controlString = NULL) {
  
  ## TODO: include proportion of variation explained
  
  dat <- ungroup(countTable)
  if(!is.null(controlString)) {
    dat <- filter(dat, grepl(controlString, gene_id, ignore.case = TRUE))
  }
  
  d <- .pcaInternal(dat, countType = "counts")
  if("norm_counts" %in% names(countTable)) {
    d <- rbind(d, .pcaInternal(dat, countType = "norm_counts"))
  }

  ggplot(d, aes(x = PC1, y = PC2)) + 
    geom_point(aes(colour = treatment, shape = replicate), size = 4) +
    facet_wrap(count_type ~ library, labeller = "label_both") +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL) +
    theme_bw() 
}


memcrispr.controlGuideCorrelation <- function(countsTable) {
  
  tmp <-
   ungroup(countsTable) %>%
   filter(grep('control', gene_id)) %>%
   mutate(GuideID = gsub(GuideID, pattern = "([0-9]+_)", replacement = "")) %>%
   select(GuideID, NormCounts, Replicate, Treatment, Library) %>%
   unite(grouping, Replicate, Treatment, Library, sep = "-") %>%
   arrange(GuideID) %>%
   spread(key = grouping, value = norm_counts)
  heatmap(cor(as.matrix(as.data.frame(tmp)[,2:9])))

  ggplot( melt(reorder_cormat(cor(as.matrix(as.data.frame(tmp)[,2:9])))) ) +
    geom_tile(aes(x = Var1, y = Var2, fill = value)) + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation")
}
  
reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }