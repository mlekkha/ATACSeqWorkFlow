# These are the R functions written by Manu lab (University of North Dakota; Biology Department) to analyze ATAC-Seq data
# These functions count Tn5 cuts per nucleotide and correct for centering/sequence bias.

# Authors: Trevor Long, Manu Manu
# Last Modified: 06/20/24

### Main function used to count Tn5 cuts at nucleotides
# samplesheet: a csv file with at least 3 columns (SampleID = sample names, BAMfile = path to sample BAM file, reads = library size)
# chrom: The chromosome of the region of interest (character; i.e "chr7")
# regionstart: The starting position of the region of interest (integer; i.e 35118303)
# regionend: The ending position of the region of interest (integer; i.e 35118803)
# biascorrect: Should dataframe include Tn5 sequebce bias-corrected counts (boolean; i.e TRUE)
# Fbiastable: Output of BagFoot MakeBiasCorrectionTableBAM() function with only + strand reads. This is only used in HINT correction
# Rbiastable: Output of BagFoot MakeBiasCorrectionTableBAM() function with only - strand reads. This is only used in HINT correction
# combinedbiastable: Output of BagFoot MakeBiasCorrectionTableBAM() function with all reads. This is used in manulong and TOBIAS correction
# genome: UCSC reference genome name. available genomes can be seen using BSgenome::available.genomes(). 
# centercuts: should cuts be shifted
# forwardShift: how many bp should + strand cuts be shifted
# reverseShift: how many bp should - strand cuts be shifted
# is...: filter parameters for scanbamparam
# nmer: size of kmer used in determining sequence bias. Should be the same length as used in bias table
# windowsize: size of smoothing window for HINT bias correction
# biascorrectionmethod: which bias correction to use ("manulong", "hint", "bagfoot")
# filtering: should scanbamparam filters be used
# p_value_threshold: adj. p_value threshold for minimum cut count needed to be corrected in manulong correction
# numsamplesinpool: used to calculate background lamda (do not need to change)

createATACDataTable <- function(samplesheet = "./",
                                chrom = NA,
                                regionstart = NA,
                                regionend = NA,
                                biascorrect = FALSE,
                                Fbiastable = NA,
                                Rbiastable = NA,
                                combinedbiastable = NA,
                                genome = "mm10",
                                centercuts = TRUE,
                                forwardShift = 4,
                                reverseShift = -4,
                                isPaired = T, 
                                isProperPair = T, 
                                isUnmappedQuery = F, 
                                isSecondaryAlignment = F, 
                                isNotPassingQualityControls = F, 
                                isDuplicate = F, 
                                isSupplementaryAlignment = F,
                                nmer = 6,
                                windowsize = 50,
                                biascorrectionmethod = "manulong",
                                filtering = TRUE,
                                p_value_threshold = 0.05,
                                genomeAssembled = TRUE,
                                genomeCircular = FALSE,
                                numSamplesinPool = 1){
  
  # load in required packages and functions
  require(reshape2, quietly = TRUE) || BiocManager::install("reshape2")
  require(ggplot2, quietly = TRUE) || BiocManager::install("ggplot2")
  require(zoo, quietly = TRUE) || BiocManager::install("zoo")
  require(gridExtra, quietly = TRUE) || BiocManager::install("gridExtra")
  require(GenomicAlignments, quietly = TRUE) || BiocManager::install("GenomicAlignments")
  require(BSgenome, quietly = TRUE) || BiocManager::install("BSgenome")
  genomeBS <- suppressMessages(grep(genome, available.genomes(), value = TRUE)[[1]])
  if (!require(genomeBS, character.only = TRUE)){
    library(genomeBS, character.only = TRUE)
  }
  require(GenomeInfoDb, quietly = TRUE) || BiocManager::install("GenomeInfoDb")

  # If biascorrect is false, function will count raw cuts
  if (biascorrect == FALSE){
    
    # call makeLocusGRange to create a region-of-interest grange
    Locus <- makeLocusGRange(chrom = chrom,
                             rstart = regionstart,
                             rend = regionend)
    
    # call countTn5eventsByRegion to count the cuts in the region-of-interest
    counts <- countTn5eventsByRegion(samplesheet = samplesheet,
                                     regions = Locus,
                                     genome = genome,
                                     centercuts = centercuts,
                                     forwardShift = forwardShift,
                                     reverseShift = reverseShift,
                                     isPaired = isPaired, 
                                     isProperPair = isProperPair, 
                                     isUnmappedQuery = isUnmappedQuery, 
                                     isSecondaryAlignment = isSecondaryAlignment, 
                                     isNotPassingQualityControls = isNotPassingQualityControls, 
                                     isDuplicate = isDuplicate, 
                                     isSupplementaryAlignment = isSupplementaryAlignment,
                                     filtering = filtering,
                                     genomeAssembled = genomeAssembled,
                                     genomeCircular = genomeCircular)
  
    # Sum the cuts since they are split by strand
    counts <- rowSums(counts[[1]], dims = 2)
    
    # Melt the data 
    counts <- melt(counts)
    
    # Name the columns
    colnames(counts) <- c("Location", "Sample", "Cuts")
    
    }
  
  
  # if biascorrect is set to true, the function will return bias-corrected cuts
  if (biascorrect == TRUE & biascorrectionmethod == "hint" | biascorrect == TRUE & biascorrectionmethod == "manu"){
    
    # call makeLocusGRange to create a region-of-interest grange
    Locus <- makeLocusGRange(chrom = chrom,
                             rstart = regionstart - windowsize,
                             rend = regionend + windowsize)
   
    # call countTn5eventsByRegion to count the cuts in the region-of-interest
    counts <- countTn5eventsByRegion(samplesheet = samplesheet,
                                     regions = Locus,
                                     genome = genome,
                                     centercuts = centercuts,
                                     forwardShift = forwardShift,
                                     reverseShift = reverseShift,
                                     isPaired = isPaired, 
                                     isProperPair = isProperPair, 
                                     isUnmappedQuery = isUnmappedQuery, 
                                     isSecondaryAlignment = isSecondaryAlignment, 
                                     isNotPassingQualityControls = isNotPassingQualityControls, 
                                     isDuplicate = isDuplicate, 
                                     isSupplementaryAlignment = isSupplementaryAlignment,
                                     filtering = filtering,
                                     genomeAssembled = genomeAssembled,
                                     genomeCircular = genomeCircular)
   
  
    # Plus strand #
    # pull the first element of the 3rd dimension from the countTn5eventsByRegion output matrix
    # first element of the 3rd dimension represents the cuts on the "+" strand
    counts_plus <- counts[[1]][,,1]
    
    # melt the matrix
    # melting the matrix creates 3 columns (Location, Sample, Cuts on "+" strand)
    counts_plus <- melt(counts_plus)
    
    # Name the columns for ease of navigation
    colnames(counts_plus) <- c("Location", "Sample", "Cuts")
    
    # call assignNMERSequence to assign the nmer sequence to each location for bias allocation
    counts_plus <- assignNMERSequence(counts_plus,
                                      nmer = nmer,
                                      regstart = regionstart - windowsize,
                                      regend = regionend + windowsize,
                                      chrom = chrom,
                                      genome = genome)
    
    # call assignBiasValue to assign the appropriate bias based on nmer sequence for each location
    counts_plus <- assignBiasValue(counts_plus,
                                   biastable = Fbiastable)
    
    # call biasCorrection to get bias corrected cuts
    counts_plus <- biasCorrection(counts_plus,
                                  method = biascorrectionmethod,
                                  windowsize = windowsize,
                                  p_value_threshold = p_value_threshold)
    
    
    # Minus strand #
    # pull the second element of the 3rd dimension from the countTn5eventsByRegion output matrix
    # second element of the 3rd dimension represents the cuts on the "-" strand
    counts_minus <- counts[[1]][,,2]
    
    # melt the matrix
    # melting the matrix creates 3 columns (Location, Sample, Cuts on "-" strand)
    counts_minus <- melt(counts_minus)
    
    # Name the columns for ease of navigation
    colnames(counts_minus) <- c("Location", "Sample", "Cuts")
    
    # call assignNMERSequence to assign the nmer sequence to each location for bias allocation
    counts_minus <- assignNMERSequence(counts_minus,
                                      nmer = nmer,
                                      regstart = regionstart - windowsize,
                                      regend = regionend + windowsize,
                                      chrom = chrom,
                                      genome = genome)
    
    # call assignBiasValue to assign the appropriate bias based on nmer sequence for each location
    counts_minus <- assignBiasValue(counts_minus,
                                   biastable = Rbiastable)
    
    # call biasCorrection to get bias corrected cuts
    counts_minus <- biasCorrection(counts_minus,
                                  method = biascorrectionmethod,
                                  windowsize = windowsize,
                                  p_value_threshold = p_value_threshold)
    
    
    # combine
    counts <- data.frame(Location = counts_plus$Location,
                         Sample = counts_plus$Sample,
                         Cuts = counts_plus$Cuts + counts_minus$Cuts,
                         correctedcuts = counts_plus$correctedcuts + counts_minus$correctedcuts)
    
    # get rid of extensions causes by smoothing
    counts <- counts[counts$Location %in% regionstart:regionend,]
  }
  
  if (biascorrect == TRUE & biascorrectionmethod == "bagfoot"){
    
    # call makeLocusGRange to create a region-of-interest grange
    Locus <- makeLocusGRange(chrom = chrom,
                             rstart = regionstart,
                             rend = regionend)
    
    # call countTn5eventsByRegion to count the cuts in the region-of-interest
    counts <- countTn5eventsByRegion(samplesheet = samplesheet,
                                     regions = Locus,
                                     genome = genome,
                                     centercuts = centercuts,
                                     forwardShift = forwardShift,
                                     reverseShift = reverseShift,
                                     isPaired = isPaired, 
                                     isProperPair = isProperPair, 
                                     isUnmappedQuery = isUnmappedQuery, 
                                     isSecondaryAlignment = isSecondaryAlignment, 
                                     isNotPassingQualityControls = isNotPassingQualityControls, 
                                     isDuplicate = isDuplicate, 
                                     isSupplementaryAlignment = isSupplementaryAlignment,
                                     filtering = filtering,
                                     genomeAssembled = genomeAssembled,
                                     genomeCircular = genomeCircular)
    
    
    # Sum the cuts since they are split by strand
    counts <- rowSums(counts[[1]], dims = 2)
    
    # Melt the data 
    counts <- melt(counts)
    
    # Name the columns
    colnames(counts) <- c("Location", "Sample", "Cuts")
    
    # call assignNMERSequence to assign the nmer sequence to each location for bias allocation
    counts <- assignNMERSequence(counts,
                                 nmer = nmer,
                                 regstart = regionstart,
                                 regend = regionend,
                                 chrom = chrom,
                                 genome = genome)
    
    # call assignBiasValue to assign the appropriate bias based on nmer sequence for each location
    counts <- assignBiasValue(counts,
                              biastable = combinedbiastable)
    
    # call biasCorrection to get bias corrected cuts
    counts <- biasCorrection(counts,
                             method = biascorrectionmethod,
                             windowsize = windowsize,
                             p_value_threshold = p_value_threshold)
    
  }
  
  if (biascorrect == TRUE & biascorrectionmethod == "manulong"){
    
    # call makeLocusGRange to create a region-of-interest grange
    Locus <- makeLocusGRange(chrom = chrom,
                             rstart = regionstart,
                             rend = regionend)
    
    # call countTn5eventsByRegion to count the cuts in the region-of-interest
    counts <- countTn5eventsByRegion(samplesheet = samplesheet,
                                     regions = Locus,
                                     genome = genome,
                                     centercuts = centercuts,
                                     forwardShift = forwardShift,
                                     reverseShift = reverseShift,
                                     isPaired = isPaired, 
                                     isProperPair = isProperPair, 
                                     isUnmappedQuery = isUnmappedQuery, 
                                     isSecondaryAlignment = isSecondaryAlignment, 
                                     isNotPassingQualityControls = isNotPassingQualityControls, 
                                     isDuplicate = isDuplicate, 
                                     isSupplementaryAlignment = isSupplementaryAlignment,
                                     filtering = filtering,
                                     genomeAssembled = genomeAssembled,
                                     genomeCircular = genomeCircular)
    
    
    # Sum the cuts since they are split by strand
    counts <- rowSums(counts[[1]], dims = 2)
    
    # Melt the data 
    counts <- melt(counts)
    
    # Name the columns
    colnames(counts) <- c("Location", "Sample", "Cuts")
    
    # Find Lamda of each sample
    counts <- CalculateLambdaperSample(counts, libsizefile=samplesheet, genome = genome,
                                       genomeAssembled = genomeAssembled,
                                       genomeCircular = genomeCircular, numSamplesinPool = numSamplesinPool)
    
    # Find P-Values of cut likelihood at each position
    counts <- CalculateAdjustedPValues(counts, genome = genome, numSamplesinPool = numSamplesinPool,
                                       genomeAssembled = genomeAssembled,
                                       genomeCircular = genomeCircular)
    
    # call assignNMERSequence to assign the nmer sequence to each location for bias allocation
    counts <- assignNMERSequence(counts,
                                 nmer = nmer,
                                 regstart = regionstart,
                                 regend = regionend,
                                 chrom = chrom,
                                 genome = genome)
    
    # call assignBiasValue to assign the appropriate bias based on nmer sequence for each location
    counts <- assignBiasValue(counts,
                              biastable = combinedbiastable)
    
    # call biasCorrection to get bias corrected cuts
    counts <- biasCorrection(counts,
                             method = biascorrectionmethod,
                             windowsize = windowsize,
                             p_value_threshold = p_value_threshold)
    
  }
  
  return(counts)
}

### Function to assign each position with its kmer bias value from matrix
assignBiasValue <- function(x,
                            biastable = NA){
  
  # read in bias table
  biastable <- read.table(biastable)  
  
  # initiate bias column
  x$bias <- 0
    
  # use for-loop to assign bias from bias table
  for (position in seq(length(x$nmer))){
    
    # addition to address no-call bases
    if (is.na(grep("N", x$nmer[position])[1]) == F){
      x$bias[position] <- 1
      print("N base found, one more more bases were assigned a correction of 1")
      next
    }
    
    x$bias[position] <- biastable$CorrectionFactor[rownames(biastable) == x$nmer[position]]
  }
  
 return(x)
}

### function to find kmer centered on each nucleotide
assignNMERSequence <- function(x,
                               nmer = 6,
                               regstart = NA,
                               regend = NA,
                               chrom = NA,
                               genome = "mm10"){
  
  # load in genome sequence needed to find/assign nmer sequence
  genomesequence <- suppressMessages(get(grep(genome, available.genomes(), value = T)[1]))
  
  # initiate a new column in the data table to place nmer sequence
  x$nmer <- "Empty"
  
  # use for-loop to assign nmer sequence to each location in x
  #x must have a column called "Location"
  for (location in seq(regstart,regend)){
    sequenceregion <- seq(location-floor((nmer-1)/2),
                          location+ceiling((nmer-1)/2))
    nmersequence <- as.character(genomesequence[[chrom]][min(sequenceregion):max(sequenceregion)])
    x$nmer[x$Location == location] <- nmersequence
  }
  
  return(x)
}

### function to correct Tn5 cuts at nucleotides based on kmer bias value
biasCorrection <- function(x,
                           method = "manulong",
                           windowsize = 1,
                           p_value_threshold = 0.05){
  
  if (method != "manu" & method != "hint" & method != "bagfoot" & method != "manulong"){
    stop(print("Requested bias correction method not available"))
  }
    
    # Manu method of bias correction divides the observed cuts at a location by the average bias within the provided window at that location
    if (method == "manu"){
      for (sampl in unique(x$Sample)){
        # create the average bias in a provided window at each location
        x$averagebias[x$Sample == sampl] <- rollmean(x$bias[x$Sample == sampl], k=windowsize, fill = 0)
        x$correctedcuts[x$Sample == sampl] <- (x$Cuts[x$Sample == sampl])/(x$averagebias[x$Sample == sampl])
      }
    }
  
    if (method == "hint"){
      for (sampl in unique(x$Sample)){
        x$rollmeancounts[x$Sample == sampl] <- rollmean(x$Cuts[x$Sample == sampl], k=windowsize, fill = 0)
        x$rollsumbias[x$Sample == sampl] <- rollsum(x$bias[x$Sample == sampl], k=windowsize, fill = 0)
        x$normalizedbias[x$Sample == sampl] <- (x$bias[x$Sample == sampl])/(x$rollsumbias[x$Sample == sampl])
        x$correctedcuts[x$Sample == sampl] <- (x$Cuts[x$Sample == sampl] +1)/(x$rollmeancounts[x$Sample == sampl] * x$normalizedbias[x$Sample == sampl] +1)
    }
  }
  
  if (method == "bagfoot"){
    x$correctedcuts <- (x$Cuts)*(x$bias)
  }
  
  if (method == "manulong"){
    
    x$passes_threshold <- FALSE
    x$correctedcuts <- x$Cuts
    log10_alpha <- log10(p_value_threshold)
    naturallog_alpha <- log10_alpha/log10(exp(1))
    
   for (location in unique(x$Location)){
     y <- x$adj_P_values[x$Location == location]
     for (qvalue in 1:length(y)){
       if (y[qvalue] <= naturallog_alpha){
         x$passes_threshold[x$Location == location] <- TRUE
         break
       }
     }
   }
    

    for (row in 1:length(x$passes_threshold)){
      if (x$passes_threshold[row] == TRUE){
        x$correctedcuts[row] <- x$Cuts[row] * x$bias[row]
      }
    }
  }
  
  return(x)
}

### Main counting workhorse function
countTn5eventsByRegion <- function(samplesheet = "./",
                                   regions = NULL, 
                                   genome = "mm10", 
                                   centercuts = FALSE,
                                   libsizefile = NULL,
                                   namerows = TRUE,
                                   forwardShift = 4,
                                   reverseShift = -4,
                                   isPaired = T, 
                                   isProperPair = T, 
                                   isUnmappedQuery = F, 
                                   isSecondaryAlignment = F, 
                                   isNotPassingQualityControls = F, 
                                   isDuplicate = F, 
                                   isSupplementaryAlignment = F,
                                   filtering = T,
                                   genomeAssembled = T,
                                   genomeCircular = F
                                  )
{
  
  ###Load libraries needed
  require(GenomicAlignments, quietly = TRUE) || BiocManager::install("GenomicAlignments")
  require(GenomicRanges, quietly = TRUE) || BiocManager::install("GenomicRanges")
  require(Rsamtools, quietly = TRUE) || BiocManager::install("Rsamtools")
  require(GenomicFeatures, quietly = TRUE) || BiocManager::install("GenomicFeatures")
  
  # check that regions is not NULL, that it is a GRanges object and that
  # it is not empty
  stopifnot(!is.null(regions), 
            class(regions) == "GRanges", 
            length(regions) > 0
           )


  # set the offset according to whether we are centering the Tn5 cuts or
  # not
  if (centercuts)
  {

    offset <- c(forwardShift,reverseShift) 

  } else {

    offset <- c(0,0) 

  }
  
  # Read in BAM files and sample info
  BAMfiles <- read.csv(samplesheet)
  BAMfilenames <- list(BAMfiles$BAMfile)
  BAMfilenames <- BAMfilenames[[1]]

  # how many samples
  numSamples <- length(BAMfilenames)

  # sample numbers
  sampleNumbers <- BAMfiles$SampleID

  # get chromosome information
  chrominfo <- suppressMessages(getBSGenomeInfo(genome = genome, genomeAssembled = genomeAssembled,
                               genomeCircular = genomeCircular))
  
  ###Subset chromosome sizes
  chromsizes <- chrominfo$seqlengths
  
  ###Subset chromosome names
  chromnames <- row.names(chrominfo)
  
  
  #We don't allow the below to run since it take too much memory and
  #crashes R. The calling function must call chromosome by chromosome.
  #If no regions were specified, then they are the chromosomes
  if (is.null(regions)) {

      regions <- GRanges(seqnames = chromnames,
                         ranges = IRanges(
                                    start = rep(1,length(chromsizes)),
                                    end = chromsizes
                                         )
                        )                 

  } 

  #Initialize Tn5counts and naming vector
  numRegions <- length(regions)

  # declare or create the list to be returned
  Tn5counts <- vector("list", 0)

  # determine the chromosome names, regions starts, and region ends
  chrs <- as.character(seqnames(regions))     
  rstarts <- start(regions)
  rends <- end(regions)

  # populate it with zero matrices of the correct dimensions
  for (i in 1:numRegions) 
  {

      chr <- chrs[i]
      rstart <- rstarts[i]
      rend <- rends[i]
      rlength <- rend - rstart + 1
      rname <- paste(chr, ":", rstart, "-", rend, sep="")

      Tn5counts[[rname]] <- array(data=0, 
                                  dim=c(rlength, numSamples, 2)
                                 ) 

      ## i.e. assign the sample names to the 2nd dimension (columns) of
      ## the 3D counts array
      dimnames(Tn5counts[[rname]])[[2]] <- sampleNumbers

      ## Name the + or - strand count columns as "+" or "-" respectively
      dimnames(Tn5counts[[rname]])[[3]] <- c("+", "-")

      ## Name the rows by genomic coordinates
      if (namerows)
          dimnames(Tn5counts[[rname]])[[1]] <- seq(rstart,rend)

  }               


  ## Create a variable to use as a filtering element wile converting BAM
  ## file to GRange if the bamflags have not already been passed as an
  ## argument
  if (filtering == TRUE){
    bamflags <- scanBamFlag(isPaired = isPaired, 
                            isProperPair = isProperPair, 
                            isUnmappedQuery = isUnmappedQuery, 
                            isSecondaryAlignment = isSecondaryAlignment, 
                            isNotPassingQualityControls = isNotPassingQualityControls, 
                            isDuplicate = isDuplicate, 
                            isSupplementaryAlignment = isSupplementaryAlignment)
    
    # Construct the set of parameters for reading the BAM files, which are
    # the SAM format FLAGs, what parameters (position and strand of the
    # read),  as well as the regions (in GRanges format) we want to read in
    bamparams_region <- ScanBamParam(flag=bamflags,
                                     simpleCigar = TRUE,
                                     what=c("qwidth", "pos", "mpos", "strand"),
                                     which=regions)
    
    # Construct the set of parameters for reading the BAM files, which are
    # the SAM format FLAGs but regions don't need to be specified since we
    # will count the reads genomewide
    bamparams_all <- ScanBamParam(flag=bamflags,
                                  simpleCigar = TRUE)
  }
 
  if (filtering == FALSE){
    bamflags <- scanBamFlag(isPaired = NA, 
                            isProperPair = NA, 
                            isUnmappedQuery = NA, 
                            isSecondaryAlignment = NA, 
                            isNotPassingQualityControls = NA, 
                            isDuplicate = NA, 
                            isSupplementaryAlignment = NA)
    
    bamparams_region <- ScanBamParam(flag=bamflags,
                                     what=c("qwidth", "pos", "mpos", "strand"),
                                     which=regions)
  }


   for (BAMfilename in BAMfilenames) {

    # keep track of sample number
    samplenumber <- BAMfiles$SampleID[BAMfiles$BAMfile == BAMfilename]

    # open the file
    theBAMfile <- BamFile(file=BAMfilename,
                          index=paste(BAMfilename,"bai",sep=".")
                         )

    # open the file
    open(theBAMfile)

    # check if the BAM file is sorted. If not, print error and exit
    bamheader <- scanBamHeader(theBAMfile)
    
    # position of the "SO:coordinate" string which indicates that the BAM
    # is sorted
    issortedstringpos <- grep("SO:coordinate", 
                              bamheader$text$`@HD`[2], 
                              fixed=TRUE
                             ) 
    # If there are no matches, then raise error and quit
    if (length(issortedstringpos) == 0)
        stop(paste("BAM file",BAMfilename,"is not sorted"))


    # Read in the reads. with_which_label instructs readGAlignments to
    # include an additional field containing the names of the regions
    # specified using the "which" parameter of  

    reads <- as.data.frame(readGAlignments(theBAMfile, 
                                           with.which_label = TRUE,
                                           param=bamparams_region
                                          )
                          )                

    readsbyregion <- split(reads, reads$which_label)
    regionnames <- names(readsbyregion)

    # Loop over the reads from the region, adjust the Tn5 cut position
    # for reads mapping to - strand, and do the counting.
    for (i in 1:length(readsbyregion))
    {

        regionreads <- readsbyregion[[i]]
        numreads <- nrow(regionreads)

        # skip to the next region if this one contains no reads
        if (numreads == 0)
            next

        # The Tn5 cut position is the "start" coordinate for reads that map
        # to the + strand
        tn5cuts <- regionreads$start

        # The Tn5 cut position is the "end" coordinate for reads that map
        # to the - strand
        maptominus <- regionreads$strand == "-"
        tn5cuts[maptominus] <- regionreads$end[maptominus]

        # Determine the strand's index in the counts, 1 for + and 2 for -
        strand <- as.numeric(maptominus) + 1

        # Adjust the cut location based on the offsets
        tn5cuts <- tn5cuts + offset[strand]

        # Determine region's coordinates to filter Tn5 cuts lying outside
        regionstring <- strsplit(regionnames[i], ":")[[1]][2]  
        regioncoords <- as.numeric(strsplit(regionstring, "-")[[1]])

        # Do the counting, making sure that the position actually lies in
        # the region
        for (j in 1:numreads)
        {

            # make sure that the cut lies within the region - scanBam
            # includes reads that overlap the region so only one cut need
            # be within the region
            if ( (tn5cuts[j] >= regioncoords[1]) &&
                    (tn5cuts[j] <= regioncoords[2]) )
                
            {

                m <- strand[j]
                l <- tn5cuts[j] - regioncoords[1] + 1 

                # counting
                Tn5counts[[regionnames[i]]][l,samplenumber,m] <-
                        Tn5counts[[regionnames[i]]][l,samplenumber,m] + 1

            }

        }

    } # from loop over regions

  close(theBAMfile)

  } # from loop over BAMfilenames
        
  return(Tn5counts)

}

### Function to take input locations and convert to a genomic range
makeLocusGRange <- function(chrom = NA,
                            rstart = NA,
                            rend = NA){
  
  # Input error handling
  if (class(chrom) != "character"){
    stop(paste("chrom must be a character (i.e 'chr7')"))
  }
  
  # load in required package to convert dataframe into grange/grangelist
  require(GenomicAlignments, quietly = TRUE) || BiocManager::install("GenomicAlignments")

  # Create a dataframe of input information
  # if vectors are given in the inputs, first row of dataframe will represent the first element of the input vectors
  Locus <- data.frame(seqnames = chrom, start = rstart, end = rend)
  
  # Convert dataframe into a Grange/grangelist
  Locus <- makeGRangesFromDataFrame(Locus)
  
  # return the grange
  return(Locus)
}

### function to calculate background signal
CalculateLambdaperSample <- function(x,
                                     libsizefile=samplesheet,
                                     genome = "mm10",
                                     genomeAssembled = TRUE,
                                     genomeCircular = FALSE,
                                     numSamplesinPool = 1){
  
  # Error handling
  if (is.null(libsizefile) == TRUE){
    stop(paste0("Please provide library sizes"))
  }
  
  # Initiate a lambda column
  x$Lambda <- 0
  
  # Read in library sizes
  librarysizes <- read.csv(libsizefile)
  
  # get chromosome information
  chrominfo.clps <- suppressMessages(getBSGenomeInfo(genome = genome, genomeAssembled = genomeAssembled,
                                    genomeCircular = genomeCircular))
  
  # Find sum of genome size
  sizeofgenome.clps <- sum(chrominfo.clps$seqlengths)
  
  # Give the ability to pool samples for lambda
  sizeofgenome.clps <- sizeofgenome.clps * numSamplesinPool
  
  # Divide library sizes by the size of the genome to find lambda of each sample
  librarysizes$Reads <- librarysizes$Reads/sizeofgenome.clps
  
  # error handling
  if (length(librarysizes$Sample) != length(unique(x$Sample))){
    stop(paste0("The number of samples in provided library sizes sheet does not match number of counted samples. Please check that sample names are correct."))
  }
  
  if (sum(unique(librarysizes$Sample) %in% unique(x$Sample)) != length(unique(x$Sample))){
    stop(paste0("Sample names in library sizes sheets do not match BAM sample names. Please check to make sure the names match."))
  }
  
  # Assign the lambda to each sample
  for (sampl in unique(librarysizes$SampleID)){
    x$Lambda[x$Sample == sampl] <- librarysizes$Reads[librarysizes$SampleID == sampl]
  }
  
  # Return to dataframe
  return(x)
}

### function to calculate adjusted p values
CalculateAdjustedPValues <- function(x,
                                     genome = "mm10",
                                     numSamplesinPool = 1,
                                     genomeAssembled = TRUE,
                                     genomeCircular = FALSE){
  
  # libraries
  require(GenomeInfoDb, quietly = TRUE) || BiocManager::install("GenomeInfoDb")
  
  # Initilize a p_value column
  x$P_Value <- 0
  
  # Assign p_values to each location based on lambda
  for (samp in unique(x$Sample)){
    x$P_Value[x$Sample == samp] <- dpois(x$Cuts[x$Sample == samp], x$Lambda[x$Sample == samp][1], log = TRUE)
  }
  
  # get chromosome info
  chrominfo.capv <- suppressMessages(getBSGenomeInfo(genome = genome, genomeAssembled = genomeAssembled,
                                    genomeCircular = genomeCircular))
  
  # Find sum of genome size
  sizeofgenome.capv <- sum(chrominfo.capv$seqlengths)
  
  # Give the ability to pool samples for correction
  sizeofgenome.capv <- sizeofgenome.capv * numSamplesinPool
  
  # add the log p values and the log size of the genome
  x$adj_P_values <- x$P_Value + log(sizeofgenome.capv)
  
  return(x)
}

### Function to get genome information
getBSGenomeInfo <- function(genome = "mm10",
                            genomeAssembled = TRUE,
                            genomeCircular = FALSE){
  
  # read in chromosome information from BSgenome
  genomeBS <- suppressMessages(grep(genome, available.genomes(), value = T)[1])
  chrominfo <- as.data.frame(SeqinfoForBSGenome(genomeBS))
  
  # subset out assembled and non-circular chromosomes if desired
  if (genomeAssembled){
    chrominfo <- chrominfo[grep("_", row.names(chrominfo), invert = T),]
  }
  chrominfo <- chrominfo[chrominfo$isCircular == genomeCircular,]
  
  # return chromosome info
  return(chrominfo)
  
}

### Function to identify extended regions of depleted signal and the associated sequences
IdentifyFootprints <- function(x,
                               genome = "mm10",
                               chrom = "chr7",
                               regionstart = 1,
                               regionend = 2,
                               minWidth = 10,
                               maxWidth = 50,
                               manualThreshold = NA){
  
  # Load in libraries needed
  require(dplyr, quietly = TRUE) || BiocManager::install("dplyr")
  require(ggplot2, quietly = TRUE) || BiocManager::install("ggplot2")
  require(BSgenome, quietly = TRUE) || BiocManager::install("BSgenome")

  # Summarize the cuts
  footprintdf <- x %>%
    group_by(Location) %>%
    summarise(total = sum(correctedcuts))
  
  # Subset the region of interest
  subfoot <- footprintdf[footprintdf$Location %in% regionstart:regionend,]
  
  
  ### THIS SECTION IS FOR AUTOMATICALLY FINDING A THRESHOLD. IF DATA HAS LOW READS, SET A MANUAL THRESHOLD ###
  
  # Find peaks that pass can be distinguished from background
  sigloc <- unique(x$Location[x$passes_threshold == T &
                                x$Location %in% regionstart:regionend])
  
  # Find the cut rate lambda in the region of interest
  # Only calculating lambda off of peaks that can be distinguished from background
  lambda <- mean(subfoot$total[subfoot$Location %in% sigloc])
  
  #############################################################################################################
  
  # manually set threshold if wanted
  if (!is.na(manualThreshold)){
    lambda <- manualThreshold
  }
  
  # Subset the locations that pass the lambda
  sigsubfoot <- subfoot[subfoot$total >= lambda,]
  
  # Initiate a footprint regions df
  loc <- data.frame(starts = 0,
                    ends = 0,
                    widths = 0)
  
  # Save the region if the space between two significant peaks are far enough apart
  for (row in 1:(length(sigsubfoot$Location) - 1)){
    n <- sigsubfoot$Location[row + 1] - sigsubfoot$Location[row]
    if (n >= minWidth){
      region <- data.frame(starts = sigsubfoot$Location[row],
                           ends = sigsubfoot$Location[row + 1],
                           widths = (sigsubfoot$Location[row + 1]) - sigsubfoot$Location[row])
      loc <- rbind(loc, region)
    }
  }
  
  # Get rid of the initiation row
  loc <- loc[-1,]
  
  # Subset out the regions that are too far apart
  loc <- loc[loc$widths <= maxWidth,]
  
  # Load in the desired genome (mm10 right now)
  genomesequence <- suppressMessages(get(grep(genome, available.genomes(), value = T)[1]))
  
  # Initiate an empty sequence column
  loc$seq <- "Empty"
  
  # Find the sequence of each region
  for (row in 1:length(loc$seq)){
    sequenceoffoot <- as.character(genomesequence[[chrom]][loc$starts[row]:loc$ends[row]])
    loc$seq[row] <- sequenceoffoot
  }
  
  # Create a binding dataframe for annotation in plotting or to just return
  color <- colors()
  numbersteps <- length(loc$starts)
  numbersteps <- (numbersteps * 10) + 8
  colorstep <- seq(9,numbersteps, by = 10)
  color <- color[colorstep]
  
  # create a dataframe with binding info
  binding <- data.frame(xmin = loc$starts,
                        xmax = loc$ends,
                        seq = loc$seq,
                        ymin=0,
                        ymax=Inf,
                        colors = color)
  
  # Plot
  bindingp <- ggplot(data = footprintdf[footprintdf$Location %in%
                                          regionstart:regionend,],
                     aes(x=Location, y=total)) +
    geom_col() + 
    annotate(geom = "rect", xmin = binding$xmin, xmax = binding$xmax,
             ymin = binding$ymin, ymax = binding$ymax,
             fill = binding$colors, alpha = 0.2) +
    theme_bw() +
    theme(axis.title=element_text(size=10, family="Helvetica")) + 
    theme(axis.text=element_text(size=8, family="Helvetica", colour="black")) + 
    theme(plot.margin = unit(c(0,0.25,0,0.125), "in"),
          plot.background=element_blank()) + 
    theme(panel.grid = element_blank()) +
    ylab("Cuts") +
    xlab("Position")
  
  
  # Return the binding data
  binding_summary <<- binding[,-c(4:6)]
  print(paste0(lambda, " Tn5 cuts set as threshold."))
  plot(bindingp)
  return(binding_summary)
  
}

### Function to match TF PWMs to sequences
MatchJASPARMotifs <- function(JASPARmatrix = "./", 
                              sequence = "ATCG",
                              min.score = "80%",
                              prior.params = c(A=0.25, C=0.25,
                                               T=0.25, G=0.25),
                              motifname = NA){
  
  # Libraries
  require(TFBSTools, quietly = TRUE) || BiocManager::install("TFBSTools")
  require(Biostrings, quietly = TRUE) || BiocManager::install("Biostrings")
  require(MatrixGenerics, quietly = TRUE) || BiocManager::install("MatrixGenerics")
  
  # Read in JASPAR motif matrix (MUST BE IN JASPAR FORMAT)
  motifs <- readJASPARMatrix(JASPARmatrix)
  
  # subset motifs if only individual motif is to be tested
  if (!is.na(motifname)){
    
    motifs <- motifs[c(grep(motifname, motifs@listData, ignore.case = T))]
    
  }
  
  # initialize final data matrix
  isMatchedfinal <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(isMatchedfinal) <- c("start", "end", "width", "seq", "score", "TF")
  
  # For-loop to match all motifs to sequence
  for ( i in 1:length(motifs) ){
    
    # check to make sure there's not a mistake in the PWM
    colhits <- max(colSums(motifs[[i]]@profileMatrix))
    misnumber <- sum(colSums(motifs[[i]]@profileMatrix) != colhits)
    if ( misnumber != 0 ) {
      
      print(paste0(motifs[[i]]@name, " did not have a rectangular PFM and will not be tested. Please check motif integrity"))
      misnumber <- 0
      next
      
    }
    
    # converts hit matrix into PWM
    thispwm <- PWM(motifs[[i]]@profileMatrix,
                   type = "log2probratio",
                   prior.params = prior.params)
    
    # Finds score of one PWM
    isMatched <- matchPWM(thispwm,
                          as.character(sequence),
                          min.score=min.score,
                          with.score=TRUE)
    
    # converts results to data frame
    isMatcheddf <- as.data.frame(isMatched)
    
    # returns message if no motifs found
    if (i == length(motifs) & length(isMatchedfinal[,1]) == 0){
      
      return(paste0("No TF motifs found on sequence with threshold of ", min.score))
      
    }
    
    # if there was no match, go to next PWM
    if (length(isMatcheddf$start) == 0){
      
      next
      
    }
    
    # Find what percentage of consensus score the PWM match is
    isMatcheddf$score <- isMatched@elementMetadata$score
    
    # Save TF name
    isMatcheddf$TF <- motifs[[i]]@name
    
    # Save PWM match to a final data frame
    isMatchedfinal <- rbind(isMatchedfinal, isMatcheddf)
    
    # convert to a list to save the tested sequence as the element name (This should make it easier to do multiple sequences)
    isMatchedlist <- list(isMatchedfinal)
    names(isMatchedlist) <- as.character(sequence)
    
  }
  
  # Return
  return(isMatchedlist)

}

# function to shift BAM files for aggregate footprinting
shiftBAMs <- function(indir = "./",
                      outdirname = "shifted/",
                      seqlev = "chr14",
                      genome = "hg38",
                      positiveShift = 4L,
                      negativeShift = -4L,
                      isPaired = TRUE, 
                      isProperPair = TRUE, 
                      isUnmappedQuery = FALSE, 
                      isSecondaryAlignment = FALSE, 
                      isNotPassingQualityControls = FALSE, 
                      isDuplicate = FALSE, 
                      isSupplementaryAlignment = FALSE){
  
  # libraries
  require(ATACseqQC, quietly = TRUE) || BiocManager::install("ATACseqQC")
  require(ChIPpeakAnno, quietly = TRUE) || BiocManager::install("ChIPpeakAnno")
  require(BSgenome, quietly = TRUE) || BiocManager::install("BSgenome")
  genomeBS <- suppressMessages(grep(genome, available.genomes(), value = TRUE)[[1]])
  if (!require(genomeBS, character.only = TRUE)){
    library(genomeBS, character.only = TRUE)
  }
  require(stringr, quietly = TRUE) || BiocManager::install("stringr")
  
  # message
  print("Beginning BAM shifting")
  
  # outpath
  dir.create(paste0(indir, outdirname))
  outpath <- paste0(indir, outdirname)
  
  # collect bams
  bamList <- normalizePath(list.files(indir,
                                      pattern = ".bam$",
                                      full.names = T))
  
  # make seqlev grange for subsetting bam
  which <- suppressMessages(get(grep(genome, available.genomes(), value = T)[1]))
  which <- as(seqinfo(which)[seqlev], "GRanges")
  
  # initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(bamList), style = 3)
  pbcount <- 0
  setTxtProgressBar(pb, pbcount)
  
  # shift bams 
  for (bfile in bamList){
    
    # iteration counter
    pbcount <- pbcount + 1
    
    # get and make file name
    namefile <- str_split(basename(bfile), "[.]")[[1]][1]
    namefile <- paste0(namefile, "_", seqlev, "_shifted.bam")
    
    # read in subset of bam
    bam <- readBamFile(bfile,
                       which = which,
                       asMates = TRUE,
                       bigFile = TRUE,
                       flag = scanBamFlag(isPaired = isPaired, 
                                          isProperPair = isProperPair, 
                                          isUnmappedQuery = isUnmappedQuery, 
                                          isSecondaryAlignment = isSecondaryAlignment, 
                                          isNotPassingQualityControls = isNotPassingQualityControls, 
                                          isDuplicate = isDuplicate, 
                                          isSupplementaryAlignment = isSupplementaryAlignment))
    
    # shift
    suppressMessages(shiftGAlignmentsList(bam,
                         outbam = paste0(outpath, namefile),
                         positive = positiveShift,
                         negative = abs(negativeShift)))
    
    # progress bar
    setTxtProgressBar(pb, pbcount)
    
  }
  
  return(paste0("Shifting done! Shifted BAM files can be found in ", outpath))
  
}
