# Functions and classes migrated from the soGGi Bioconductor package (deprecated).
# Original authors: Gopuraja Dharmalingam, Tom Carroll.
# Only the components required by profileplyr are included.

# ChIPprofile S4 class -------------------------------------------------------

#' ChIPprofile S4 class
#'
#' An S4 class representing ChIP-seq profiles over genomic intervals.
#' Extends \code{RangedSummarizedExperiment} with a \code{params} slot
#' storing the parameters used to generate the profiles.
#'
#' @slot params A list of parameters used when generating the profile.
#' @exportClass ChIPprofile
setClass("ChIPprofile", contains = "RangedSummarizedExperiment",
         slots = c(params = "list"))

#' Combine ChIPprofile objects
#'
#' @param x A ChIPprofile object
#' @param ... Additional ChIPprofile objects
#' @return A ChIPprofile object
#' @export
setMethod("c", "ChIPprofile",
          function(x, ...) {
            assayList <- lapply(list(x, ...), function(x) assays(x)[[1]])
            subsetProfile <- SummarizedExperiment(assayList, rowRanges = rowRanges(x))
            metadata(subsetProfile)$names <- unlist(lapply(list(x, ...), function(x) metadata(x)$name))
            metadata(subsetProfile)$AlignedReadsInBam <- unlist(lapply(list(x, ...), function(x) metadata(x)$AlignedReadsInBam))
            return(new("ChIPprofile", subsetProfile, params = x@params))
          })


# regionPlot -----------------------------------------------------------------

#' @importFrom GenomicAlignments readGAlignments qwidth
#' @importFrom Biostrings matchPWM reverseComplement
#' @importFrom Rsamtools BamFile index indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom rtracklayer import.bw
#' @importFrom IRanges runValue shiftApply
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BiocGenerics "seqlengths<-"
#' @importFrom S4Vectors Rle
#' @importFrom stats spline
regionPlot <- function(bamFile, testRanges, samplename = NULL, nOfWindows = 100,
                       FragmentLength = 150, style = "point", distanceAround = NULL,
                       distanceUp = NULL, distanceDown = NULL,
                       distanceInRegionStart = NULL, distanceOutRegionStart = NULL,
                       distanceInRegionEnd = NULL, distanceOutRegionEnd = NULL,
                       paired = FALSE, normalize = "RPM", plotBy = "coverage",
                       removeDup = FALSE, verbose = TRUE, format = "bam",
                       seqlengths = NULL, forceFragment = NULL, method = "bin",
                       genome = NULL, cutoff = 80, downSample = NULL,
                       minFragmentLength = NULL, maxFragmentLength = NULL) {
  if (!verbose) {
    suppressMessages(runRegionPlot())
  }
  result <- runRegionPlot(bamFile, testRanges, samplename, nOfWindows, FragmentLength,
                          style, distanceAround, distanceUp, distanceDown,
                          distanceInRegionStart, distanceOutRegionStart,
                          distanceInRegionEnd, distanceOutRegionEnd, paired, normalize,
                          plotBy, removeDup, format, seqlengths, forceFragment, method,
                          genome, cutoff, downSample, minFragmentLength, maxFragmentLength)
  return(result)
}

runRegionPlot <- function(bamFile, testRanges, samplename = NULL, nOfWindows = 100,
                          FragmentLength = 150, style = "point", distanceAround = NULL,
                          distanceUp = NULL, distanceDown = NULL,
                          distanceInRegionStart = NULL, distanceOutRegionStart = NULL,
                          distanceInRegionEnd = NULL, distanceOutRegionEnd = NULL,
                          paired = FALSE, normalize = "RPM", plotBy = "coverage",
                          removeDup = FALSE, format = "bam", seqlengths = NULL,
                          forceFragment = NULL, method = "bin", genome = NULL,
                          cutoff = 80, downSample = NULL, minFragmentLength = NULL,
                          maxFragmentLength = NULL) {

  if (format == "bam") {
    if (file.exists(bamFile) & length(index(BamFile(bamFile))) == 0) {
      message("Creating index for ", bamFile)
      indexBam(bamFile)
      message("..done")
    }
  }

  testRanges <- GetGRanges(testRanges)
  ## Check parameters
  if (style != "percentOfRegion") {
    if (is.null(distanceAround)) {
      distanceAround = 1500
    }
    if (is.null(distanceUp)) {
      distanceUp <- distanceAround
    }
    if (is.null(distanceDown)) {
      distanceDown <- distanceAround
    }
  } else {
    if (is.null(distanceAround)) {
      distanceAround = 100
    }
    if (is.null(distanceUp)) {
      distanceUp <- distanceAround
    }
    if (is.null(distanceDown)) {
      distanceDown <- distanceAround
    }
  }

  if (is.null(distanceInRegionStart)) {
    distanceInRegionStart = 750
  }
  if (is.null(distanceOutRegionStart)) {
    distanceOutRegionStart = 1500
  }
  if (is.null(distanceInRegionEnd)) {
    distanceInRegionStart = 750
  }
  if (is.null(distanceOutRegionEnd)) {
    distanceInRegionStart = 1500
  }

  ## Initialize empty matrices and parameters for collecting coverage analysis
  ## Find maximum distance to use for filtering out of bounds extended GRanges
  if (style == "region" | style == "regionandpoint") {
    posRegionStartMat <- NULL
    posRegionEndMat <- NULL
    negRegionStartMat <- NULL
    negRegionEndMat <- NULL
    RegionsMat <- NULL
    maxDistance <- max(distanceOutRegionStart, distanceOutRegionEnd)
    distanceUpStart <- distanceOutRegionStart
    distanceDownEnd <- distanceOutRegionEnd
  }

  if (style == "point") {
    PosRegionMat <- NULL
    NegRegionMat <- NULL
    RegionsMat <- NULL
    whatIsMax <- max(distanceAround, distanceUp, distanceDown)
    maxDistance <- whatIsMax
    distanceUpStart <- distanceUp
    distanceDownEnd <- distanceDown
  }

  if (style == "percentOfRegion") {
    maxDistance <- round((distanceAround / 100) * width(testRanges))
    RegionsMat <- NULL
    distanceUpStart <- NULL
    distanceDownEnd <- NULL
  }
  totalReads <- NA

  if (format == "bam") {
    message("Reading Bam header information...", appendLF = FALSE)
    allchrs <- names(scanBamHeader(bamFile)[[1]]$targets)
    lengths <- as.vector(scanBamHeader(bamFile)[[1]]$targets)
    names(lengths) <- allchrs
    message("..Done")
  }

  if (format == "bigwig") {
    message("Importing BigWig...", appendLF = FALSE)
    genomeCov <- import.bw(bamFile, as = "RleList")
    if (is.null(seqlengths)) {
      seqlengths(genomeCov) <- unlist(lapply(genomeCov, length))
    } else {
      seqlengths(genomeCov)[match(names(lengths), names(genomeCov))] <- lengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }

  if (format == "pwm") {
    bamFile <- pwmToCoverage(bamFile, genome, min = cutoff, removeRand = FALSE)
    format <- "rlelist"
  }

  if (format == "granges") {
    genomeCov <- coverage(bamFile)
    format <- "rlelist"
  }

  if (format == "rlelist") {
    message("Importing rlelist", appendLF = FALSE)
    genomeCov <- bamFile
    if (is.null(seqlengths)) {
      seqlengths(genomeCov) <- unlist(lapply(genomeCov, length))
    } else {
      seqlengths(genomeCov)[match(names(lengths), names(genomeCov))] <- lengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }

  if (style != "percentOfRegion") {
    message("Filtering regions which extend outside of genome boundaries...", appendLF = FALSE)
    testRangeNames <- unique(seqnames(testRanges))
    temptestranges <- GRanges()
    for (i in 1:length(testRangeNames)) {
      perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
      temptestranges <- c(temptestranges, perchrRanges[end(perchrRanges) + maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                       & start(perchrRanges) - maxDistance > 0])
    }
  }
  if (style == "percentOfRegion") {
    message("Filtering regions which extend outside of genome boundaries...", appendLF = FALSE)
    testRangeNames <- unique(seqnames(testRanges))
    temptestranges <- GRanges()
    for (i in 1:length(testRangeNames)) {
      perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
      perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
      temptestranges <- c(temptestranges, perchrRanges[end(perchrRanges) + perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                       & start(perchrRanges) - perChrMaxDistance > 0])
      perChrMaxDistance <- perChrMaxDistance[end(perchrRanges) + perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                             & start(perchrRanges) - perChrMaxDistance > 0]
      distanceUpStart <- c(distanceUpStart, perChrMaxDistance)
    }
    distanceDownEnd <- distanceUpStart
  }
  message("..Done")
  message("Filtered ", length(testRanges) - length(temptestranges), " of ", length(testRanges), " regions")
  testRanges <- temptestranges
  temptestranges <- NULL

  message("Splitting regions by Watson and Crick strand..", appendLF = FALSE)
  mcols(testRanges) <- cbind(mcols(testRanges), data.frame(giID = paste0("giID", seq(1, length(testRanges)))))
  strand(testRanges[strand(testRanges) == "*"]) <- "+"
  testRangesPos <- testRanges[strand(testRanges) == "+"]
  testRangesNeg <- testRanges[strand(testRanges) == "-"]
  message("..Done")
  if (style == "percentOfRegion") {
    distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
    distanceDownEndPos <- distanceUpStartPos
    distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
    distanceDownEndNeg <- distanceUpStartNeg
    message("..Done")
  } else {
    distanceUpStartPos <- distanceUpStart
    distanceDownEndPos <- distanceDownEnd
    distanceUpStartNeg <- distanceUpStart
    distanceDownEndNeg <- distanceDownEnd
    message("..Done")
  }

  if (style == "region") {
    message("Filtering regions which are smaller than windows into region...", appendLF = FALSE)
    testRangesPos <- testRangesPos[(end(testRangesPos) - distanceInRegionEnd) - (start(testRangesPos) + distanceInRegionStart) > nOfWindows]
    testRangesNeg <- testRangesNeg[(end(testRangesNeg) - distanceInRegionStart) - (start(testRangesNeg) + distanceInRegionEnd) > nOfWindows]
    message("..Done")
  }

  message("Found ", length(testRangesPos), " Watson strand regions")
  message("Found ", length(testRangesNeg), " Crick strand regions")

  message("Extending regions..", appendLF = FALSE)
  exttestRanges <- c(GRanges(seqnames(testRangesPos), IRanges(start(testRangesPos) - distanceUpStartPos, end(testRangesPos) + distanceDownEndPos)),
                     GRanges(seqnames(testRangesNeg), IRanges(start(testRangesNeg) - distanceDownEndNeg, end(testRangesNeg) + distanceUpStartNeg)))
  message("...done")

  reducedExtTestRanges <- reduce(exttestRanges)

  if (!removeDup) {
    Param <- ScanBamParam(which = GRanges(seqnames = seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]), IRanges(start = start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]), end = end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
  } else {
    Param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), which = GRanges(seqnames = seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]), IRanges(start = start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]), end = end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
  }

  if (format == "bam") {
    message("Reading tags from ", bamFile, appendLF = FALSE)
    totalReads <- 10^6

    if (paired == FALSE) {
      total <- readGAlignments(bamFile, param = Param)
      message("..Done.\nRead in ", length(total), " reads")

      if (is.null(FragmentLength)) {
        FragmentLength <- getShifts(total, lengths, shiftWindowStart = 1, shiftWindowEnd = 400)
      }

      message("Extending reads to fragmentlength of ", FragmentLength, appendLF = FALSE)
      temp <- resize(as(total, "GRanges"), FragmentLength, "start")
      message("..done")
    }

    if (paired == TRUE) {
      gaPaired <- readGAlignments(bamFile,
                                  param = ScanBamParam(what = c("mpos"),
                                                       flag = scanBamFlag(isProperPair = TRUE, isFirstMateRead = TRUE),
                                                       mapqFilter = 30))
      tempPos <- GRanges(seqnames(gaPaired[strand(gaPaired) == "+"]),
                         IRanges(
                           start = start(gaPaired[strand(gaPaired) == "+"]),
                           end = mcols(gaPaired[strand(gaPaired) == "+"])$mpos
                           + qwidth(gaPaired[strand(gaPaired) == "+"])))
      tempNeg <- GRanges(seqnames(gaPaired[strand(gaPaired) == "-"]),
                         IRanges(
                           start = mcols(gaPaired[strand(gaPaired) == "-"])$mpos,
                           end = end(gaPaired[strand(gaPaired) == "-"])))
      temp <- c(tempPos, tempNeg)
      message("..Done.\nRead in ", length(temp), " reads")
      if (removeDup) {
        message("Removing duplicates")
        beforeDupR <- length(temp)
        temp <- unique(temp)
        AfterDupR <- length(temp)
        message("Removed ", beforeDupR - AfterDupR, " duplicates")
      }

      if (!is.null(minFragmentLength)) {
        temp <- temp[width(temp) > minFragmentLength]
      }
      if (!is.null(maxFragmentLength)) {
        temp <- temp[width(temp) < maxFragmentLength]
      }
      if (!is.null(forceFragment)) {
        message("Forcing fragments to be centred and set to ", forceFragment, "..", appendLF = FALSE)
        temp <- resize(temp, forceFragment, "center")
        message("..done")
      }
      message("..done")
    }

    message("Calculating coverage..", appendLF = FALSE)
    seqlengths(temp)[match(names(lengths), names(seqlengths(temp)))] <- lengths
    if (!is.null(downSample)) {
      if (downSample < 1 & downSample > 0) {
        temp <- temp[sample(length(temp), round(length(temp)) * downSample), ]
      } else if (downSample > 1) {
        temp <- temp[sample(length(temp), downSample), ]
      }
    }

    genomeCov <- coverage(temp)
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..done")
  }
  chromosomes <- seqlevels(genomeCov)

  if (style == "point") {
    testRangesPos <- resize(testRangesPos, 1, "center")
    testRangesNeg <- resize(testRangesNeg, 1, "center")
    RangesPos <- GRanges(seqnames(testRangesPos), IRanges(start(testRangesPos) - distanceUpStart, start(testRangesPos) + distanceDownEnd), strand = Rle("+", length(testRangesPos)), mcols(testRangesPos))
    RangesNeg <- GRanges(seqnames(testRangesNeg), IRanges(end(testRangesNeg) - distanceDownEnd, end(testRangesNeg) + distanceUpStart), strand = Rle("-", length(testRangesNeg)), mcols(testRangesNeg))
    message("Calculating coverage across regions\nCalculating per contig. ")

    for (c in 1:length(chromosomes)) {
      message(paste0("contig: ", c))
      if (length(RangesPos[seqnames(RangesPos) %in% chromosomes[c]]) > 0) {
        PosRegionMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesPos[seqnames(RangesPos) %in% chromosomes[c]])]), ncol = mean(width(RangesPos)), byrow = TRUE)
        rownames(PosRegionMat) <- RangesPos[seqnames(RangesPos) %in% chromosomes[c]]$giID
      }
      if (length(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]) > 0) {
        NegRegionMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]])])), ncol = mean(width(RangesNeg)), byrow = TRUE)
        rownames(NegRegionMat) <- RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]$giID
      }
      RegionsMat <- rbind(RegionsMat, PosRegionMat, NegRegionMat)
    }
    message("Creating ChIPprofile.")

    profileMat <- RegionsMat
    colnames(profileMat) <- c(paste0("Point_Centre", seq(0 - distanceUpStart, -1)), "Point_Centre", paste0("Point_Centre", seq(1, distanceDownEnd)))
    filteredRanges <- c(RangesPos, RangesNeg)
    profileSample <- SummarizedExperiment(profileMat, rowRanges = filteredRanges[match(rownames(profileMat), filteredRanges$giID)])

    if (is.null(samplename)) {
      if (format %in% c("rlelist", "pwm", "granges")) {
        metadata(profileSample) <- list(names = c("Sample"))
      } else {
        metadata(profileSample) <- list(names = c(bamFile), AlignedReadsInBam = totalReads)
      }
    } else {
      metadata(profileSample) <- list(names = samplename, AlignedReadsInBam = totalReads)
    }

    paramList <- list("nOfWindows" = nOfWindows,
                      "style" = style,
                      "samplename" = samplename,
                      "nOfWindows" = nOfWindows,
                      "FragmentLength" = FragmentLength,
                      "distanceAround" = distanceAround,
                      "distanceUp" = distanceUp,
                      "distanceDown" = distanceDown,
                      "distanceInRegionStart" = distanceInRegionStart,
                      "distanceInRegionEnd" = distanceInRegionEnd,
                      "distanceOutRegionStart" = distanceOutRegionStart,
                      "distanceOutRegionEnd" = distanceOutRegionEnd,
                      "paired" = paired,
                      "normalize" = normalize,
                      "plotBy" = plotBy,
                      "removeDup" = removeDup,
                      "format" = format,
                      "seqlengths" = seqlengths,
                      "forceFragment" = forceFragment,
                      "method" = method,
                      "genome" = genome,
                      "cutoff" = cutoff,
                      "minFragmentLength" = minFragmentLength,
                      "maxFragmentLength" = maxFragmentLength,
                      "downSample" = downSample)
    return(new("ChIPprofile", profileSample, params = paramList))
  }

  if (style == "percentOfRegion") {

    if (method == "spline") {
      grWidths <- width(testRangesPos)
      Flanks <- round(grWidths * ((distanceAround) / 100))
      RangesPos <- GRanges(seqnames(testRangesPos), IRanges(start(testRangesPos) - Flanks, end(testRangesPos) + Flanks), strand = Rle("+", length(testRangesPos)), mcols(testRangesPos))
      grWidths <- width(testRangesNeg)
      Flanks <- round(grWidths * ((distanceAround) / 100))
      RangesNeg <- GRanges(seqnames(testRangesNeg), IRanges(start(testRangesNeg) - Flanks, end(testRangesNeg) + Flanks), strand = Rle("+", length(testRangesNeg)), mcols(testRangesNeg))

      matPos <- NULL
      matNeg <- NULL
      testRangesPosNew <- GRanges()
      testRangesNegNew <- GRanges()
      message(paste0("Calculating splines for regions.\nProcessing per contig"))

      for (c in 1:length(chromosomes)) {
        message(paste0("contig: ", i))
        if (any(seqnames(RangesPos) == chromosomes[c])) {
          testRangesPosNew <- c(testRangesPosNew, RangesPos[seqnames(RangesPos) == chromosomes[c]])
          matPos = c(matPos, list(t(viewApply(
            Views(genomeCov[names(genomeCov) == chromosomes[c]][[1]],
                  ranges(RangesPos[seqnames(RangesPos) == chromosomes[c]])),
            function(x) spline(x, n = (2 * (nOfWindows * ((distanceAround) / 100))) + nOfWindows)$y))))
        }

        if (any(seqnames(RangesNeg) == chromosomes[c])) {
          testRangesNegNew <- c(testRangesNegNew, RangesNeg[seqnames(RangesNeg) == chromosomes[c]])
          matNeg = c(matNeg, list(t(viewApply(
            Views(genomeCov[names(genomeCov) == chromosomes[c]][[1]],
                  ranges(RangesNeg[seqnames(RangesNeg) == chromosomes[c]])),
            function(x) spline(x, n = (2 * (nOfWindows * ((distanceAround) / 100))) + nOfWindows)$y))[, ((2 * (nOfWindows * ((distanceAround) / 100))) + nOfWindows):1]))
        }
      }

      message("Creating ChIPprofile")
      if (!is.null(matPos)) {
        matPos <- do.call(rbind, matPos)
      }
      if (!is.null(matNeg)) {
        matNeg <- do.call(rbind, matNeg)
      }
      meansMat <- rbind(matPos, matNeg)
      allRanges <- c(testRangesPosNew, testRangesNegNew)
      rownames(meansMat) <- allRanges$giID
      profileMat <- meansMat[order(rownames(meansMat)), ]
      colnames(profileMat) <- c(paste0("Start-", seq(1, (nOfWindows * ((distanceAround) / 100)))),
                                paste0("Start+", seq(1, nOfWindows)),
                                paste0("End+", seq(1, (nOfWindows * ((distanceAround) / 100)))))

      profileSample <- SummarizedExperiment(profileMat, rowRanges = allRanges[match(rownames(profileMat), allRanges$giID)])

      if (is.null(samplename)) {
        if (format %in% c("rlelist", "pwm", "granges")) {
          metadata(profileSample) <- list(names = c("Sample"))
        } else {
          metadata(profileSample) <- list(names = c(bamFile), AlignedReadsInBam = totalReads)
        }
      } else {
        metadata(profileSample) <- list(names = samplename, AlignedReadsInBam = totalReads)
      }

      paramList <- list("nOfWindows" = nOfWindows,
                        "style" = style,
                        "samplename" = samplename,
                        "nOfWindows" = nOfWindows,
                        "FragmentLength" = FragmentLength,
                        "distanceAround" = distanceAround,
                        "distanceUp" = distanceUp,
                        "distanceDown" = distanceDown,
                        "distanceInRegionStart" = distanceInRegionStart,
                        "distanceInRegionEnd" = distanceInRegionEnd,
                        "distanceOutRegionStart" = distanceOutRegionStart,
                        "distanceOutRegionEnd" = distanceOutRegionEnd,
                        "paired" = paired,
                        "normalize" = normalize,
                        "plotBy" = plotBy,
                        "removeDup" = removeDup,
                        "format" = format,
                        "seqlengths" = seqlengths,
                        "forceFragment" = forceFragment,
                        "method" = method,
                        "genome" = genome,
                        "cutoff" = cutoff,
                        "minFragmentLength" = minFragmentLength,
                        "maxFragmentLength" = maxFragmentLength,
                        "downSample" = downSample)
      return(new("ChIPprofile", profileSample, params = paramList))
    }

    if (method == "bin") {
      meansListNeg <- vector("numeric")
      meansListPos <- vector("numeric")
      grListWindowsPos <- GRanges()
      grListWindowsNeg <- GRanges()

      message("Making windows.")

      if (length(testRangesPos) > 0) {
        grWidths <- width(testRangesPos)
        windows <- floor(grWidths %/% nOfWindows)
        extraForWindows <- grWidths %% nOfWindows
        extraForFlankWindows <- grWidths %% (nOfWindows * ((distanceAround) / 100))
        addToWindow <- 0
        startPos <- start(testRangesPos) - distanceUpStartPos
        rem <- rep(0, length(extraForFlankWindows))
        rem2 <- NULL

        message("Windowing positive 5' flanking ")
        for (i in 1:(nOfWindows * ((distanceAround) / 100))) {
          rem2 <- rem + ((extraForFlankWindows >= i) + 0)
          grListWindowsPos <- c(grListWindowsPos, GRanges(seqnames(testRangesPos), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesPos$giID))
          rem <- rem2
        }

        startPos <- start(testRangesPos)
        rem <- rep(0, length(extraForWindows))
        rem2 <- NULL
        message("Windowing positive regions ")
        for (i in 1:(nOfWindows)) {
          rem2 <- rem + ((extraForWindows >= i) + 0)
          grListWindowsPos <- c(grListWindowsPos, GRanges(seqnames(testRangesPos), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesPos$giID))
          rem <- rem2
        }

        rem <- rep(0, length(extraForFlankWindows))
        rem2 <- NULL
        startPos <- end(testRangesPos)
        message("Windowing positive 3' flank ")
        for (i in 1:(nOfWindows * ((distanceAround) / 100))) {
          rem2 <- rem + ((extraForFlankWindows >= i) + 0)
          grListWindowsPos <- c(grListWindowsPos, GRanges(seqnames(testRangesPos), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesPos$giID))
          rem <- rem2
        }

        grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
      }

      if (length(testRangesNeg) > 0) {
        grWidths <- width(testRangesNeg)
        windows <- floor(grWidths %/% nOfWindows)
        extraForWindows <- grWidths %% nOfWindows
        extraForFlankWindows <- grWidths %% (nOfWindows * ((distanceAround) / 100))
        addToWindow <- 0
        startPos <- start(testRangesNeg) - distanceDownEndNeg
        rem <- rep(0, length(extraForFlankWindows))
        rem2 <- NULL
        message("Windowing negative 5' flank ")
        for (i in 1:(nOfWindows * ((distanceAround) / 100))) {
          rem2 <- rem + ((extraForFlankWindows >= i) + 0)
          grListWindowsNeg <- c(grListWindowsNeg, GRanges(seqnames(testRangesNeg), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesNeg$giID))
          rem <- rem2
        }

        startPos <- start(testRangesNeg)
        rem <- rep(0, length(extraForWindows))
        rem2 <- NULL
        message("Windowing negative regions ")
        for (i in 1:(nOfWindows)) {
          rem2 <- rem + ((extraForWindows >= i) + 0)
          grListWindowsNeg <- c(grListWindowsNeg, GRanges(seqnames(testRangesNeg), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesNeg$giID))
          rem <- rem2
        }
        rem <- rep(0, length(extraForFlankWindows))
        rem2 <- NULL
        startPos <- end(testRangesNeg)
        message("Windowing negative 3' flank ")
        for (i in 1:(nOfWindows * ((distanceAround) / 100))) {
          rem2 <- rem + ((extraForFlankWindows >= i) + 0)
          grListWindowsNeg <- c(grListWindowsNeg, GRanges(seqnames(testRangesNeg), IRanges(
            (startPos) + rem + (windows * (i - 1)),
            startPos + (windows * i) - 1 + rem2), giID = testRangesNeg$giID))
          rem <- rem2
        }
        grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
      }
      grListWindows <- list(grListWindowsPos, grListWindowsNeg)
      message("..done\n")

      message(paste0("Calculating bin scores for regions.\nProcessing per contig"))
      for (c in 1:length(chromosomes)) {
        message(paste0("contig: ", c))
        message("Processing inner region windows in ", chromosomes[c])
        covPerPeakPos <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]], ranges(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]))
        doubleTempPos <- viewMeans(covPerPeakPos)
        names(doubleTempPos) <- as.vector(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]$giID)
        meansListPos <- c(meansListPos, doubleTempPos)
        covPerPeakNeg <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]], ranges(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]))
        doubleTempNeg <- viewMeans(covPerPeakNeg)
        names(doubleTempNeg) <- as.vector(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]$giID)
        meansListNeg <- c(meansListNeg, doubleTempNeg)
        message("..done")
        message("Processing flanking windows in ", chromosomes[c])

        tempstartRegionRangesPosMat <- NULL
        tempendRegionRangesPosMat <- NULL
        tempstartRegionRangesNegMat <- NULL
        tempendRegionRangesNegMat <- NULL
      }

      meansPos <- matrix(meansListPos,
                         ncol = ((nOfWindows * ((distanceAround) / 100)) * 2) + nOfWindows,
                         byrow = TRUE)
      if (nrow(meansPos) > 0) {
        rownames(meansPos) <- matrix(names(meansListPos), ncol = ((nOfWindows * ((distanceAround) / 100)) * 2) + nOfWindows,
                                     byrow = TRUE)[, 1]
      }
      meansNeg <- matrix(meansListNeg,
                         ncol = ((nOfWindows * ((distanceAround) / 100)) * 2) + nOfWindows,
                         byrow = TRUE)[, (((nOfWindows * ((distanceAround) / 100)) * 2) + nOfWindows):1]
      if (nrow(meansNeg) > 0) {
        rownames(meansNeg) <- matrix(names(meansListNeg), ncol = ((nOfWindows * ((distanceAround) / 100)) * 2) + nOfWindows,
                                     byrow = TRUE)[, 1]
      }
      meansMat <- rbind(meansPos, meansNeg)
      profileMat <- meansMat[order(rownames(meansMat)), ]
    }

    colnames(profileMat) <- c(paste0("Start-", seq(1, (nOfWindows * ((distanceAround) / 100)))),
                              paste0("Start+", seq(1, nOfWindows)),
                              paste0("End+", seq(1, (nOfWindows * ((distanceAround) / 100)))))
    filteredRanges <- c(testRangesPos, testRangesNeg)
    profileSample <- SummarizedExperiment(profileMat, rowRanges = filteredRanges[match(rownames(profileMat), filteredRanges$giID)])

    if (is.null(samplename)) {
      if (format %in% c("rlelist", "pwm", "granges")) {
        metadata(profileSample) <- list(names = c("Sample"))
      } else {
        metadata(profileSample) <- list(names = c(bamFile), AlignedReadsInBam = totalReads)
      }
    } else {
      metadata(profileSample) <- list(names = samplename, AlignedReadsInBam = totalReads)
    }

    paramList <- list("nOfWindows" = nOfWindows,
                      "style" = style,
                      "samplename" = samplename,
                      "nOfWindows" = nOfWindows,
                      "FragmentLength" = FragmentLength,
                      "distanceAround" = distanceAround,
                      "distanceUp" = distanceUp,
                      "distanceDown" = distanceDown,
                      "distanceInRegionStart" = distanceInRegionStart,
                      "distanceInRegionEnd" = distanceInRegionEnd,
                      "distanceOutRegionStart" = distanceOutRegionStart,
                      "distanceOutRegionEnd" = distanceOutRegionEnd,
                      "paired" = paired,
                      "normalize" = normalize,
                      "plotBy" = plotBy,
                      "removeDup" = removeDup,
                      "format" = format,
                      "seqlengths" = seqlengths,
                      "forceFragment" = forceFragment,
                      "method" = method,
                      "genome" = genome,
                      "cutoff" = cutoff,
                      "minFragmentLength" = minFragmentLength,
                      "maxFragmentLength" = maxFragmentLength,
                      "downSample" = downSample)
    return(new("ChIPprofile", profileSample, params = paramList))
  }

  if (style == "region") {

    message("Defining flanks of regions..", appendLF = FALSE)
    startRegionRangesPos <- GRanges(seqnames(testRangesPos), IRanges(start(testRangesPos) - distanceOutRegionStart, start(testRangesPos) + distanceInRegionStart), strand = Rle("+", length(testRangesPos)), mcols(testRangesPos))
    endRegionRangesPos <- GRanges(seqnames(testRangesPos), IRanges(end(testRangesPos) - distanceInRegionEnd, end(testRangesPos) + distanceOutRegionEnd), strand = Rle("+", length(testRangesPos)), mcols(testRangesPos))
    startRegionRangesNeg <- GRanges(seqnames(testRangesNeg), IRanges(end(testRangesNeg) - distanceInRegionStart, end(testRangesNeg) + distanceOutRegionStart), strand = Rle("+", length(testRangesNeg)), mcols(testRangesNeg))
    endRegionRangesNeg <- GRanges(seqnames(testRangesNeg), IRanges(start(testRangesNeg) - distanceOutRegionEnd, start(testRangesNeg) + distanceInRegionEnd), strand = Rle("+", length(testRangesNeg)), mcols(testRangesNeg))

    testRangesPos <- GRanges(seqnames(testRangesPos), IRanges(start(testRangesPos) + distanceInRegionStart, end(testRangesPos) - distanceInRegionEnd), strand = Rle("+", length(testRangesPos)), mcols(testRangesPos))
    testRangesNeg <- GRanges(seqnames(testRangesNeg), IRanges(start(testRangesNeg) + distanceInRegionEnd, end(testRangesNeg) - distanceInRegionStart), strand = Rle("+", length(testRangesNeg)), mcols(testRangesNeg))
    message("...Done")

    meansList <- vector("numeric")
    grListWindowsPos <- GRanges()
    grListWindowsNeg <- GRanges()
    message("Making windows")

    if (length(testRangesPos) > 0) {
      grWidths <- width(testRangesPos)
      windows <- floor(grWidths %/% nOfWindows)
      extraForWindows <- grWidths %% nOfWindows
      addToWindow <- 0
      startPos <- start(testRangesPos)
      rem <- rep(0, length(extraForWindows))
      rem2 <- NULL
      message("Windowing positive regions ")
      for (i in 1:(nOfWindows)) {
        rem2 <- rem + ((extraForWindows >= i) + 0)
        grListWindowsPos <- c(grListWindowsPos, GRanges(seqnames(testRangesPos), IRanges(
          (startPos) + rem + (windows * (i - 1)),
          startPos + (windows * i) - 1 + rem2), giID = testRangesPos$giID))
        rem <- rem2
      }
      grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]

      if (length(testRangesNeg) > 0) {
        grWidths <- width(testRangesNeg)
        windows <- floor(grWidths %/% nOfWindows)
        extraForWindows <- grWidths %% nOfWindows
        addToWindow <- 0
        startNeg <- start(testRangesNeg)
        rem <- rep(0, length(extraForWindows))
        rem2 <- NULL
        message("Windowing negative regions ")
        for (i in 1:(nOfWindows)) {
          rem2 <- rem + ((extraForWindows >= i) + 0)
          grListWindowsNeg <- c(grListWindowsNeg, GRanges(seqnames(testRangesNeg), IRanges(
            (startNeg) + rem + (windows * (i - 1)),
            startNeg + (windows * i) - 1 + rem2), giID = testRangesNeg$giID))
          rem <- rem2
        }
        grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
      }
      grListWindows <- c(grListWindowsPos, grListWindowsNeg)

      message(paste0("Calculating bin scores for regions and per base pair for flanks.\nProcessing per contig"))
      for (c in 1:length(chromosomes)) {
        message(paste0("contig: ", c))
        covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]], ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
        doubleTemp <- viewMeans(covPerPeak)
        names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$giID)
        meansList <- c(meansList, doubleTemp)

        tempstartRegionRangesPosMat <- NULL
        tempendRegionRangesPosMat <- NULL
        tempstartRegionRangesNegMat <- NULL
        tempendRegionRangesNegMat <- NULL

        if (length(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]) > 0) {
          tempstartRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]])]), ncol = mean(width(startRegionRangesPos)), byrow = TRUE)
          rownames(tempstartRegionRangesPosMat) <- startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]$giID
        }
        if (length(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]) > 0) {
          tempendRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]])]), ncol = mean(width(endRegionRangesPos)), byrow = TRUE)
          rownames(tempendRegionRangesPosMat) <- endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]$giID
        }
        if (length(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]) > 0) {
          tempstartRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]])])), ncol = mean(width(startRegionRangesNeg)), byrow = TRUE)
          rownames(tempstartRegionRangesNegMat) <- rev(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]$giID)
        }
        if (length(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]) > 0) {
          tempendRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]])])), ncol = mean(width(endRegionRangesNeg)), byrow = TRUE)
          rownames(tempendRegionRangesNegMat) <- rev(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]$giID)
        }

        posRegionStartMat <- rbind(posRegionStartMat, tempstartRegionRangesPosMat)
        posRegionEndMat <- rbind(posRegionEndMat, tempendRegionRangesPosMat)
        negRegionStartMat <- rbind(negRegionStartMat, tempstartRegionRangesNegMat)
        negRegionEndMat <- rbind(negRegionEndMat, tempendRegionRangesNegMat)
        tempstartRegionRangesPosMat <- NULL
        tempendRegionRangesPosMat <- NULL
        tempstartRegionRangesNegMat <- NULL
        tempendRegionRangesNegMat <- NULL
        message("..done")
      }

      message("Creating ChIPprofile")
      AllRegionStart <- rbind(posRegionStartMat, negRegionStartMat)
      AllRegionEnd <- rbind(posRegionEndMat, negRegionEndMat)
      meansMat <- matrix(meansList, ncol = nOfWindows, byrow = TRUE)
      rownames(meansMat) <- matrix(names(meansList), ncol = nOfWindows, byrow = TRUE)[, 1]
      start <- cbind(seq(1, length(colMeans(AllRegionStart))), colMeans(AllRegionStart))
      mid <- cbind(max(start[, 1]) + seq(1, length(colMeans(meansMat))) * nOfWindows, colMeans(meansMat))
      end <- cbind(max(mid[, 1]) + seq(1, length(colMeans(AllRegionEnd))), colMeans(AllRegionEnd))
      profileMat <- cbind(AllRegionStart[order(rownames(AllRegionStart)), ],
                          meansMat[order(rownames(meansMat)), ],
                          AllRegionEnd[order(rownames(AllRegionEnd)), ])
      colnames(profileMat) <- c(paste0("Region_Start", seq(0 - distanceOutRegionStart, -1)), "Region_Start", paste0("Region_Start", seq(1, distanceInRegionStart)),
                                paste0(seq(1, nOfWindows), "%_ofRegion"),
                                paste0("Region_End", seq(0 - distanceInRegionEnd, -1)), "Region_End", paste0("Region_End", seq(1, distanceOutRegionEnd)))
      filteredRanges <- c(testRangesPos, testRangesNeg)
      profileSample <- SummarizedExperiment(profileMat, rowRanges = filteredRanges[match(rownames(profileMat), filteredRanges$giID)])
      print(format)

      if (is.null(samplename)) {
        if (format %in% c("rlelist", "pwm", "granges")) {
          metadata(profileSample) <- list(names = c("Sample"))
        } else {
          metadata(profileSample) <- list(names = c(bamFile), AlignedReadsInBam = totalReads)
        }
      } else {
        metadata(profileSample) <- list(names = samplename, AlignedReadsInBam = totalReads)
      }

      paramList <- list("nOfWindows" = nOfWindows,
                        "style" = style,
                        "samplename" = samplename,
                        "nOfWindows" = nOfWindows,
                        "FragmentLength" = FragmentLength,
                        "distanceAround" = distanceAround,
                        "distanceUp" = distanceUp,
                        "distanceDown" = distanceDown,
                        "distanceInRegionStart" = distanceInRegionStart,
                        "distanceInRegionEnd" = distanceInRegionEnd,
                        "distanceOutRegionStart" = distanceOutRegionStart,
                        "distanceOutRegionEnd" = distanceOutRegionEnd,
                        "paired" = paired,
                        "normalize" = normalize,
                        "plotBy" = plotBy,
                        "removeDup" = removeDup,
                        "format" = format,
                        "seqlengths" = seqlengths,
                        "forceFragment" = forceFragment,
                        "method" = method,
                        "genome" = genome,
                        "cutoff" = cutoff,
                        "minFragmentLength" = minFragmentLength,
                        "maxFragmentLength" = maxFragmentLength,
                        "downSample" = downSample)
      return(new("ChIPprofile", profileSample, params = paramList))
    }
  }
}


# Internal helpers for runRegionPlot -----------------------------------------

GetGRanges <- function(LoadFile, AllChr = NULL, ChrOfInterest = NULL, simple = FALSE,
                       sepr = "\t", simplify = FALSE) {
  if (is(LoadFile, "GRanges")) {
    RegionRanges <- LoadFile
    if (simplify) {
      RegionRanges <- GRanges(seqnames(RegionRanges), ranges(RegionRanges))
    }
  } else {
    if (is(LoadFile, "character")) {
      RangesTable <- read.delim(LoadFile, sep = sepr, header = TRUE, comment.char = "#")
    } else if (is(LoadFile, "matrix")) {
      RangesTable <- as.data.frame(LoadFile)
    } else {
      RangesTable <- as.data.frame(LoadFile)
    }
    Chromosomes <- as.vector(RangesTable[, 1])
    Start <- as.numeric(as.vector(RangesTable[, 2]))
    End <- as.numeric(as.vector(RangesTable[, 3]))
    RegionRanges <- GRanges(seqnames = Chromosomes, ranges = IRanges(start = Start, end = End))
    if (simple == FALSE) {
      if (ncol(RangesTable) > 4) {
        ID <- as.vector(RangesTable[, 4])
        Score <- as.vector(RangesTable[, 5])
        if (ncol(RangesTable) > 6) {
          Strand <- rep("*", nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[, -c(1:6)])
          mcols(RegionRanges) <- cbind(ID, Score, Strand, RemainderColumn)
        } else {
          mcols(RegionRanges) <- cbind(ID, Score)
        }
      }
    }
  }
  if (!is.null(AllChr)) {
    RegionRanges <- RegionRanges[seqnames(RegionRanges) %in% AllChr]
    seqlevels(RegionRanges, pruning.mode = "coarse") <- AllChr
  }
  if (!is.null(ChrOfInterest)) {
    RegionRanges <- RegionRanges[seqnames(RegionRanges) == ChrOfInterest]
  }
  return(RegionRanges)
}

RleSumAny <- function(e1, e2) {
  if (!requireNamespace("chipseq", quietly = TRUE))
    stop("Package 'chipseq' is required for automatic fragment length detection. Install with BiocManager::install('chipseq')")
  len <- length(e1)
  stopifnot(len == length(e2))
  x1 <- runValue(e1); s1 <- cumsum(runLength(e1))
  x2 <- runValue(e2); s2 <- cumsum(runLength(e2))
  .Call("rle_sum_any",
        as.integer(x1), as.integer(s1),
        as.integer(x2), as.integer(s2),
        as.integer(len),
        PACKAGE = "chipseq")
}

runGetShifts <- function(reads, ChrLengths, ChrOfInterestshift,
                         shiftWindowStart = 1, shiftWindowEnd = 400) {
  reads <- reads
  ChrLengths <- seqlengths(reads)
  PosCoverage <- coverage(IRanges(start(reads[strand(reads) == "+"]), start(reads[strand(reads) == "+"])),
                          width = ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
  NegCoverage <- coverage(IRanges(end(reads[strand(reads) == "-"]), end(reads[strand(reads) == "-"])),
                          width = ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
  message("Calculating shift for ", ChrOfInterestshift, "\n")
  ShiftsTemp <- shiftApply(seq(shiftWindowStart, shiftWindowEnd), PosCoverage, NegCoverage, RleSumAny, verbose = TRUE)
  return(ShiftsTemp)
}

getShifts <- function(reads, ChrLengths, shiftWindowStart = 1, shiftWindowEnd = 400) {
  if (is.character(reads)) {
    reads <- readGAlignments(reads)
  }
  shiftMat <- do.call(cbind, bplapply(names(ChrLengths), function(x)
    runGetShifts(reads[seqnames(reads) %in% x], ChrLengths, x,
                 shiftWindowStart = 1, shiftWindowEnd = 400)))
  cc_scores <- (rowSums(shiftMat)[1] - rowSums(shiftMat)) / rowSums(shiftMat)[1]
  return(cc_scores)
}


# PWM helpers -----------------------------------------------------------------

#' PWM hits as an RleList
#'
#' Converts a PWM matrix to an RleList of motif hit coverage, suitable for
#' use as input to \code{regionPlot()}.
#'
#' @param pwm A PWM matrix object.
#' @param genome A BSgenome object.
#' @param min PWM score cutoff as a percentage of maximum score (default "70\%").
#' @param removeRand Logical; remove contigs with "random" in their name.
#' @param chrsOfInterest Character vector of chromosomes to include.
#' @return An RleList of motif density per base pair.
#' @export
pwmToCoverage <- function(pwm, genome, min = "70%", removeRand = FALSE,
                          chrsOfInterest = NULL) {
  allchrs <- seqnames(genome)
  if (!is.null(allchrs)) {
    allchrs <- allchrs[allchrs %in% chrsOfInterest]
  }
  if (removeRand) {
    allchrs <- allchrs[!grepl("random", allchrs, ignore.case = TRUE)]
  }
  intergerList <- lapply(allchrs, function(x) pwmHitAsCoverage(pwm, genome, min, x))
  myrle <- RleList(intergerList, compress = FALSE)
  names(myrle) <- allchrs
  myrle
}

pwmHitAsCoverage <- function(pwm, genome, min, chrofinterest) {
  posMotifs <- matchPWM(pwm, genome[[chrofinterest]], min.score = min)
  negMotifs <- matchPWM(reverseComplement(pwm), genome[[chrofinterest]], min.score = min)
  if (length(posMotifs) > 0) {
    rleMotifHitPos <- coverage(GRanges(seqnames = chrofinterest, ranges(posMotifs), strand = "+"),
                               width = length(genome[[chrofinterest]]))
  } else {
    rleMotifHitPos <- RleList(rep(0, length(genome[[chrofinterest]])))
  }
  if (length(negMotifs) > 0) {
    rleMotifHitNeg <- coverage(GRanges(seqnames = chrofinterest, ranges(negMotifs), strand = "-"),
                               width = length(genome[[chrofinterest]]))
  } else {
    rleMotifHitNeg <- RleList(rep(0, length(genome[[chrofinterest]])))
  }
  rleTotal <- rleMotifHitPos + rleMotifHitNeg
  return(rleTotal[[1]])
}


# Dataset documentation ------------------------------------------------------

#' Example ChIPprofile object (Ikaros profiles)
#'
#' ChIP-seq signal profiles over peak regions from two Ikaros antibodies,
#' stored as a ChIPprofile object.
#'
#' @docType data
#' @name ik_Profiles
#' @usage data(ik_Profiles)
#' @return A ChIPprofile object
NULL

#' Example ChIPprofile object
#'
#' ChIP-seq signal profiles over gene bodies, stored as a ChIPprofile object.
#'
#' @docType data
#' @name chipExampleBig
#' @usage data(chipExampleBig)
#' @return A ChIPprofile object
NULL
