
#' Normalize individual guides within a sample to account for GC content
#' 
#' Normalization
#' 
#' @return data.table of counts
#' @importFrom stringr str_count
#' @export
normalizeWithinSamples <- function(countTable) {
  
  countTable <- 
    countTable %>%
    mutate(gc = stringr::str_count(string = seq, pattern = "G|C")) %>%
    group_by(sample_id, treatment, replicate, library, gc) %>% 
    mutate(gc_counts = as.numeric(median(counts))) %>% 
    group_by(sample_id, treatment, replicate, library) %>% 
    mutate(gc_counts = counts / (gc_counts / median(unique(gc_counts)))) %>%
    ungroup()
  
  return(countTable)
}

#' Normalizing between samples to account for varying total read counts
#' 
#' Normalization
#' 
#' @return data.table of counts
#' @export
normalizeBetweenSamples.old <- function(countTable) {
  
  ## if we have correct for GC, use those values
  if('GC_Counts' %in% colnames(countTable)) {
    countTable <- 
      group_by(countTable, sample_id, treatment, replicate, library) %>% 
      mutate(norm_counts = 100 * GC_Counts / median(GC_Counts[grep("control", gene_id)]))
  } else { ## here we use the raw counts
    countTable <- 
      group_by(countTable, sample_id, treatment, replicate, library) %>% 
      mutate(norm_counts = 100 * counts / median(counts[grep("control", gene_id, ignore.case = TRUE)]))
  }
  
  return(countTable)
}

#' Normalizing between samples to account for varying total read counts
#' 
#' Normalization
#' 
#' @return data.table of counts
#' @export
normalizeBetweenSamples.readdepth <- function(countTable) {
  
    countTable <- 
      group_by(countTable, sample_id, treatment, replicate, library) %>% 
      mutate(norm_counts = 1e7 * counts / sum(counts))

  return(countTable)
}

#' Normalizing between samples to account for varying total read counts
#' 
#' Normalization using 'median ratio method' proposed in DESeq 
#' 
#' @return data.table of counts
#' @export
memcrispr.normalizeBetweenSamples <- function(countTable) {
  
    tmpTab <- 
      group_by(countTable, gene_id) %>% 
      mutate(geom_mean = gm_mean(counts)) %>%
      ungroup() %>% 
      group_by(sample_id, treatment, replicate, library) %>%
      mutate(size_factor = median(counts / geom_mean)) %>%
      ungroup() %>%
      mutate(norm_counts = counts / size_factor) %>%
      select(-size_factor, -geom_mean)

  return(tmpTab)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}