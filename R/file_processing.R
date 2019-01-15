
readCounts <- function(file, path = NULL, sampleID = NULL, treatment = NULL, replicate = NULL) {
    countsTable <- read.table(file.path(path, file), sep="\t", stringsAsFactors=FALSE)
    ## remove the last line, this holds the counts of unmapped reads
    countsTable <- countsTable[-nrow(countsTable), ]
    return(data_frame("sample_id" = sampleID, "treatment" = treatment, "replicate" = replicate, "guide_id" = countsTable[,1], "counts" = countsTable[,3]))
}


readGuideSequences <- function(file, path = NULL, libraryName = NULL) {
    guides <- read.table(file.path(path, file), sep="\t", stringsAsFactors=FALSE)
    guides <- data_frame("Library" = libraryName, "guide_id" = .fixGuideNames(guides[,1]), "Gene" = .extractGeneSymbol(guides[,1]), "GuideSeq" = guides[,2])
    guides <- mutate(guides, GC = as.integer(lapply(gregexpr("C|G", GuideSeq), length)))
    return(guides)
}

## updated to work with the newer format and updated ids
readGuideSequences2 <- function(file, path = NULL, libraryName = NULL) {
  guides <- read.delim(file.path(path, file), sep=",", stringsAsFactors=FALSE)
  guides <- data_frame("Library" = .libFromGuideName(guides[,2]), "guide_id" = guides[,2], "Gene" = guides[,1], "GuideSeq" = guides[,3])
  guides <- mutate(guides, GC = as.integer(lapply(gregexpr("C|G", GuideSeq), length)))
  return(guides)
}

print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

readGuideSequencesSimple <- function(file) {
  guides <- read.delim(file.path(file), sep="\t", stringsAsFactors=FALSE)
  guides <- as_data_frame(guides)
  if(anyDuplicated(guides$guide_id)) {
    idx <- sort(c(which(duplicated(tmp$guide_id)), 
                  which(duplicated(tmp$guide_id, fromLast = TRUE))))
    stop("Duplicate guide_id entries found in ", file,
         "\nDuplicate entries are:\n", 
         print_and_capture(guides[idx,]),
         call. = FALSE)
  }
  return(guides)
}

.libFromGuideName <- function(guides) {
    lib <- gsub(guides, pattern = "_[0-9]+", replacement = "")
    return(lib)
}

## some guides are missing the 'rv)' in their name, which makes matching fail.
## this finds those instances and fixes them.
.fixGuideNames <- function(names) {
    names <- gsub(names, pattern = "(\\()([0-9]+)$", replacement = "\\1\\2 rv)")
    #names <- gsub(names, pattern = "([0-9]+_)", replacement = "")
    return(names)
}

.extractGeneSymbol <- function(names) {
  genes <- sub(x = names, pattern = "[0-9]*_([a-zA-Z0-9_\\.-]+)\\|?.*", replacement = "\\1")
  return(genes)
}


.readSampleSheet <- function(file, path) {
  ss <- read.table(file.path(path, file), header = TRUE)
  if(!all(c("File", "Sample", "Treatment", "Replicate") %in% colnames(ss))) { 
    stop("Expected to find the following columns in the sample sheet:\n
         File, Sample, Treatment")
  }
  return(ss)
}


#' Read text files of counts
#' 
#' This is the main data import function for the \code{MEMcrispR} package.  It 
#' expects to be provided with a directory location containing a sample sheet
#' and text files holding the observed counts for each guide, with one file per
#' sequencing library.
#' 
#' @param path File path to the folder containing the data.
#' @param sampleSheet Name of the sample sheet file. If left as \code{NULL}
#' this will default to `\code{SampleSheet.txt}'.
#' 
#' @return A data_frame containing counts for all guides across all samples.
#' @export
memcrispr.readCounts <- function(path, sampleSheet = NULL, guideLibraries = NULL) {
  
  if(is.null(sampleSheet))
    sampleSheet <- "SampleSheet.txt"
  
  ss <- .readSampleSheet(file = sampleSheet, path = path)
  
  ## read each of the counts files
  message("Reading counts")
  allCounts <- list()
  for(i in 1:nrow(ss)) {

    allCounts[[i]] <- readCounts(file = ss[i,"File"], path = path, 
                         sampleID = ss[i,"Sample"], 
                         treatment = ss[i,"Treatment"], 
                         replicate = ss[i,"Replicate"] )
  }
  allCounts <- do.call("rbind", allCounts)
  
  ## load the guide sequences
  ## we do this to include the sequences for the guide
  message("Checking sequences")
  if(is.null(guideLibraries)) {
    data('human.gecko.v2.libA', 'human.gecko.v2.libB', 
         'human.TKO.v1.base', 'human.TKO.v1.supp', 
         'human.CRiNCL.v1',
        package = 'MEMcrispR',
        envir = environment())
    lib_gRNAs <- bind_rows(human.gecko.v2.libA, human.gecko.v2.libB, 
                       human.TKO.v1.base, human.TKO.v1.supp,
                       human.CRiNCL.v1)
  } else {
    guideLibraries <- guideLibraries[file.exists(guideLibraries)]
    if(!length(guideLibraries))
      stop('Cannot find guide files')
    lib_gRNAs <- as_data_frame(
      dplyr::bind_rows(lapply(guideLibraries, read.table, 
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)))
  }
  
  finalTable <- left_join(allCounts, lib_gRNAs, by = "guide_id") %>% 
    as_data_frame() %>%
    arrange(gene_id) %>%
    mutate(sample_id = as.character(sample_id), replicate = as.integer(replicate), 
           treatment = as.integer(treatment), gene_id = as.character(gene_id), 
           library = as.character(library), guide_id = as.character(guide_id)) 
  
  return(finalTable)
}