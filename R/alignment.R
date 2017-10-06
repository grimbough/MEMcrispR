#' @importFrom Biostrings DNAStringSet BStringSet 
#' @importFrom ShortRead ShortRead writeFasta
.createReferenceFasta <- function(guideID = NULL, guideSequences = NULL, guideLibrary = NULL) {
  
  seqs <- DNAStringSet(paste0("NNNNCTTGTGGAAAGGACGAAACACCG",
                              as.character(guideSequences), 
                              "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCAC"))
  
  id <- BStringSet(gsub(pattern = " ", replacement = "*", x = paste(as.character(guideLibrary), as.character(guideID), sep = ":")))
  
  fasta <- ShortRead(id = id, sread = seqs)
  tmpFile <- tempfile()
  writeFasta(object = fasta, file = tmpFile)
  
  return(tmpFile)
  
}



.selectReference <- function(guideLibraries) {
  
  referenceTab <- list()
  if("GeCKOv2_A" %in% guideLibraries) 
    referenceTab[["GeCKOv2_A"]] <- human.gecko.v2.libA
  if("GeCKOv2_B" %in% guideLibraries) 
    referenceTab[["GeCKOv2_B"]] <- human.gecko.v2.libB    
  if("TKOv1_base" %in% guideLibraries) 
    referenceTab[["TKOv1_base"]] <- human.TKO.v1.base
  if("TKOv1_supp" %in% guideLibraries) 
    referenceTab[["TKOv1_supp"]] <- human.TKO.v1.supp  
  
  ## if we didn't find any matches, see if this is a file
  if(length(referenceTab) == 0) {
    for(f in guideLibraries) {
      if(file.exists(f)) {
        tmp <- readGuideSequencesSimple(f)
        referenceTab[[ tmp$library[1] ]] <- tmp
      }
    }
  }
  
  referenceTab <- do.call('rbind', referenceTab)
  
  return(referenceTab)
}

#' Align reads to reference
#' 
#' This function expects to be provided a folder containing FASTQ files 
#' produced by sequencing a genetic screening experiment.  It will align
#' the reads against a reference sequence generated from specified guide
#' libraries.
#' 
#' @param path Full file path to the folder containing the data.
#' @param sampleSheet Name of the file specifying the structure of the
#' experiment.  If left NULL a default of 'SampleSheet.txt' is used.
#' @param outputDir Full path to the location where output files should be 
#' written.  If left NULL the current working directory is used.
#' @param guideLibraries Names of libraries to align against. These can be
#' paths to files for custom libraries not included in the package.
#' @param ncores The number of CPU cores to be used by the Rsubread 
#' aligner. Defaults to 4.
#' 
#' @importFrom Rsubread align buildindex
#' @importFrom tools file_path_sans_ext
#' @export
metacrispr.align <- function(path, sampleSheet = NULL, outputDir = NULL, 
                          guideLibraries = c("GeCKOv2_A", "GeCKOv2_B"),
                          ncores = 4){
  
  if(is.null(sampleSheet))
    sampleSheet <- "SampleSheet.txt"
  
  ss <- .readSampleSheet(file = sampleSheet, path = path)
  
  referenceTab <- .selectReference(guideLibraries)
  reference <- .createReferenceFasta(referenceTab[['guide_id']],
                                     referenceTab[['seq']],
                                     referenceTab[['library']])
  
  idx <- tempfile()

  buildindex(basename = idx, reference = reference)
  
  files <- paste(path, ss[,'File'], sep = .Platform$file.sep)

  if(!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  bamFiles <- paste0(outputDir, .Platform$file.sep, 
                     basename(file_path_sans_ext(files, compression = TRUE)), '.bam')

  align(index=idx,
        readfile1 = files, 
        output_file = bamFiles,
        maxMismatches = 2, indels=0, 
        nsubreads = 8, type = 'dna', nthreads = ncores,
        unique = FALSE, nBestLocations = 2)
  
  # write out a modified sample sheets
  ss_bam <- ss %>% 
    mutate(File = stringr::str_replace(string = File, 
                                       pattern = ".fq.gz$|.fastq.gz$|.fq$|.fastq$", 
                                       replacement = ".bam"))
  write.table(ss_bam, file = file.path(outputDir, "SampleSheet.txt"), row.names = FALSE, quote = FALSE)
}

#' Count occurrences of each guide in aligned data.
#' 
#' This function expects to be provided a folder containing BAM files 
#' produced either by \code{\link{metacrispr.align}} or another aligner. It will
#' count the number of reads mapping to each guide in the reference and 
#' produce a table of the same structure as produced by \code{samtools idxstat}.
#' 
#' @param path Full file path to the folder containing the data.
#' @param sampleSheet Name of the file specifying the structure of the
#' experiment.  If left NULL a default of 'SampleSheet.txt' is used.
#' @param outputDir Full path to the location where output files should be 
#' written.  If left NULL the current working directory is used.
#' @param guideLibraries Names of libraries to align against.  These can be
#' paths to files for custom libraries not included in the package.
#' @param ncores The number of CPU cores to be used by the Rsubread 
#' aligner. Defaults to 4.
#' 
#' @importFrom Rsubread featureCounts 
#' @importFrom tidyr separate gather
#' @importFrom tibble as_data_frame data_frame
#' @export
metacrispr.count <- function(path, sampleSheet = NULL, outputDir = NULL, 
                          guideLibraries = c("GeCKOv2_A", "GeCKOv2_B"),
                          ncores = 4) {
  
  if(is.null(sampleSheet))
    sampleSheet <- "SampleSheet.txt"
  
  ss <- .readSampleSheet(file = sampleSheet, path = path)
  
  bamFiles <- file.path(path, ss[,'File'])
  
  referenceTable <- .selectReference(guideLibraries)
  
  
  id <- paste(as.character(referenceTable[['library']]), as.character(referenceTable[['guide_id']]), sep = ":")
  annotation <- data_frame(GeneID = as.character(id),
                           Chr = as.character(id),
                           Start = 1,
                           End = 112,
                           Strand = '*')
  
  count_stats <- featureCounts(files = bamFiles, 
                          annot.ext = annotation, 
                          nthreads = ncores,
                          countMultiMappingReads = TRUE)
  
  counts <- count_stats$counts
  
  colnames(counts) <- basename(bamFiles)
  counts2 <- cbind(guide_id = rownames(counts), as_data_frame(counts)) %>%
    as_data_frame() %>%
    separate(guide_id, into = c("library", 'guide_id'), sep = ':')
  
  # here we try to work out which reference a sample actually aligned to
  tmp <- gather(counts2, key = bam_file, value = counts, 3:ncol(counts2)) %>%
    group_by(bam_file, library) %>%
    summarize(counts = sum(counts))
  
  tmp2 <- .identifyLibFromSample(tmp)
  
  ## TODO: create output folder if it doesn't exist
  if(!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  .writeCounts(outputDir = outputDir, counts = counts2, librarySelection = tmp2)
  
  # write out a modified sample sheets
  ss_counts <- ss %>% 
    mutate(File = stringr::str_replace(string = File, pattern = ".bam", replacement = ".txt"))

  write.table(ss_counts, file = file.path(outputDir, "SampleSheet.txt"), 
              row.names = FALSE, quote = FALSE)
  
  message("Read count files written to ", outputDir)
}


.writeCounts <- function(outputDir, counts, librarySelection) {

  for(i in 3:ncol(counts)) {
    libs <-  filter(librarySelection, bam_file == colnames(counts)[i])[['library']]
    subset <- select(counts, c(1,2,i)) %>%
      filter(library %in% libs)
    tab <- data_frame(subset[['guide_id']], 112, subset[[3]], 0)
    ## TODO: include unaligned counts for QC stats
    tab <- rbind(tab, c('*', 0, 0, 0))
    
    outfile <- file.path(outputDir, paste0(tools::file_path_sans_ext(names(counts)[i]), ".txt"))
    
    write.table(tab, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t',
                file = outfile)
  }
  
}

#' function to try and identify which library(ies) were actually used for the 
#' sample, based on the number of reads that aligned.
#' @noRd
.identifyLibFromSample <- function(x) {
  
  x2 <- group_by(x, bam_file) %>%
    mutate(prop = counts / sum(counts)) %>% 
    filter(prop > 0.2) %>%
    select(bam_file, library) %>%
    ungroup()
  
}
