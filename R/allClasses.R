if(getRversion() >= "2.15.1")  utils::globalVariables(c("Sample","Signal","myTempGR","ranges_combined",
                                                        "info","range_combined","rowDataToOrderBy","cutree","anno_summary",
                                                        "kmeans","Groups","."))


#' @rdname profileplyr
#' @import GenomicRanges
#' @export
setClass("profileplyr",
         contains = "RangedSummarizedExperiment",
         slots=c(params="list",
                 sampleData="DataFrame",
                 sampleParams="DataFrame"
         ))





### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###
.profileplyr.charbound <-
  function(idx, txt, fmt)
  {
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
      msg <- paste(selectSome(orig[bad]), collapse=" ")
      stop(sprintf(fmt, msg))
    }
    idx
  }

.subsetprofileplyr <- function(x, i, j, k, ..., drop = FALSE)
{
  # if (1L != length(drop) || (!missing(drop) && drop))
  #   warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
  
  if (missing(i) && missing(j) && missing(k))
    return(x)
  
  if (!missing(k)) {
    message(k)
  }
  if (!missing(i)) {
    if (is.character(i)) {
      fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
      i <- .profileplyr.charbound(i, rownames(x), fmt)
    }
    ii <- as.vector(i)
    ans_elementMetadata <- x@elementMetadata[i, , drop=FALSE]
    ans_rowRanges <- x@rowRanges[i]
    ans_assays <- x@assays[ii,]$data
    # x <- BiocGenerics:::replaceSlots(x, ...,
    #                                    elementMetadata=ans_elementMetadata,
    #                                    rowRanges=ans_rowRanges,
    #                                    assays=ans_assays,
    #                                    check=FALSE)
    tempDou <- SummarizedExperiment(ans_assays,
                                    rowRanges=ans_rowRanges)
    
    metadata(tempDou)$info <- metadata(x)$info
    
    metadata(tempDou)$info$group_boundaries <- c(which(!duplicated(mcols(tempDou)$dpGroup))-1,length(mcols(tempDou)))
    metadata(tempDou)$sampleData <- metadata(x)$sampleData
    metadata(tempDou)$rowGroupsInUse <- metadata(x)$rowGroupsInUse
    
    x <- new("profileplyr", tempDou,params=x@params,
             sampleParams=x@sampleParams,
             sampleData=x@sampleData)
    
  
    
  }
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- .profileplyr.charbound(j, colnames(x), fmt)
    }
    ans_colData <- x@colData[j, , drop=FALSE]
    jj <- as.vector(j)
    ans_assays <- x@assays[,jj]$data
    # x <- BiocGenerics:::replaceSlots(x, ...,
    #                                    colData=ans_colData,
    #                                    assays=ans_assays,
    #                                    check=FALSE)
    
    tempDou <- SummarizedExperiment(ans_assays,
                                    rowRanges=rowRanges(x))
    
    metadata(tempDou)$info <- metadata(x)$info
    
    metadata(tempDou)$info$group_boundaries <- c(which(!duplicated(mcols(tempDou)$dpGroup))-1,length(mcols(tempDou)))
    metadata(tempDou)$sampleData <- metadata(x)$sampleData
    metadata(tempDou)$rowGroupsInUse <- metadata(x)$rowGroupsInUse
    
    x <- new("profileplyr", tempDou,params=x@params,
             sampleParams=x@sampleParams,
             sampleData=x@sampleData)
    
  }
  if (!missing(k)) {
    if (is.character(k)) {
      fmt <- paste0("<", class(x), ">[,k] index out of bounds: %s")
      k <- .profileplyr.charbound(k, assayNames(x), fmt)
    }
    metaData <- x@metadata
    metaData$sampleData <- x@metadata$sampleData[k, , drop=FALSE]
    kk <- as.vector(k)
    ans_assays <- x@assays$data[kk]
    
    # names(ans_assays) <- names(x@assays)
    # x <- BiocGenerics:::replaceSlots(x, ...,
    #                                  metadata=metaData,
    #                                  assays=ans_assays,
    #                                  check=FALSE)
    tempDou <- SummarizedExperiment(ans_assays,
                                    rowRanges=rowRanges(x))
    
    metadata(tempDou) <- metaData
    
    x <- new("profileplyr", tempDou,params=x@params,
             sampleParams=metaData$sampleData,
             sampleData=metaData$sampleData)
  }
  

  return(x)
}

#' Retrieve and set sample data in profileplyr object
#' @rdname sampleData
#' @param object A profileplyr object
#' @return  A DataFrame containing sample data
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' sampleData(object)
#' sampleData(object)$scale <- c(1,10,1)
#' @export
setGeneric("sampleData", function(object="profileplyr") standardGeneric("sampleData"))

#' @rdname sampleData
#' @return  A DataFrame containing sample data to replace current sample data
#' @export
setMethod("sampleData", "profileplyr",
          function (object)
          {
            return(object@sampleData)
          }
)


#' @rdname sampleData
setGeneric("sampleData<-", function(object, value)
  standardGeneric("sampleData<-"))

#' @param value DataFrame of sample information
#' @rdname sampleData
#' @export
setReplaceMethod("sampleData", c("profileplyr", "DataFrame"),
                 function(object, value) {
                   if(nrow(sampleData(object)) != nrow(value)) stop("Replacement sampleData must have same number of rows as current sampleData")
                   if(is.null(rownames(sampleData(object)))) stop("Replacement sampleData must rownames")
                   metadata(object)$sampleData <- value
                   object@sampleData <- value
                   names(assays(object)) <- rownames(value)
                   return(object)
                 })


sampleParams <- function (x)
          {
            return(x@sampleParams)
          }


#' Join, subset and manipulate ChIPprofile objects
#' @rdname profileplyr
#' @return  A profileplyr object
#' @export
setMethod("c", "profileplyr",
          function (x,...)
          {
            if(all(unlist(lapply(list(x,...),
                                 function(l) all(rowRanges(l)==rowRanges(x))
                                 )
                          )
                   )
               ){
            assayList <- unlist(as(lapply(list(x,...),function(l)assays(l)),"SimpleList"),use.names = FALSE)
            sampleParamsList <- lapply(list(x,...),function(l)sampleParams(l))
            sampleDataList <- lapply(list(x,...),function(l)sampleData(l))
            metaX <- metadata(x)
            newSampleParams <- do.call(rbind,sampleParamsList)
            newSampleData <- do.call(rbind,sampleDataList)
            names(assayList) <- rownames(newSampleData)
            dpAssays <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))

            x <- new("profileplyr", dpAssays,
                     params=x@params,
                     sampleParams=newSampleParams,
                     sampleData=newSampleData)
            metadata(x) <- metaX
            metadata(x)$sampleData <- newSampleData
            return(x)
            }else{
              stop("All profileplyr objects must have same rowData")
              
            }
          }
)


#' @param x profileplyr object
#' @param i An integer or character scalar indicating ranges of profileplyr object to return
#' @param j An integer or character scalar indicating columns of profileplyr object to return or a An integer or character scalar indicating which profileplyr object samples to return
#' @param k An integer or character scalar indicating samples of profileplyr object to return.
#' @param ... Additional arguments.
#' @param drop A logical whether to drop empty samples
#' @rdname profileplyr
setMethod("[", c("profileplyr", "ANY", "ANY", "ANY"),
          .subsetprofileplyr)

#' @rdname profileplyr
#' @export
setMethod("[[", c("profileplyr", "ANY", "missing"),
          function(x, i, j, ...)
          {
            subsetProfile <- SummarizedExperiment(assays(x)[[i, ...]],rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names[i]
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam[i]
            metadata(subsetProfile)$info <- c(info)
            metadata(subsetProfile)$sampleData <- sampleData[i,,drop=FALSE]
            metadata(subsetProfile)$info$group_boundaries <- c(which(!duplicated(myTempGR$dpGroup))-1,length(myTempGR))
            tempDou <- new("profileplyr", subsetProfile)
            
            return(tempDou)                        
          })

.DollarNames.profileplyr <- function(object, pattern = "")
  grep(pattern, rownames(metadata(object)$sampleData), value=TRUE)


#' Dataframe of top differentially expressed genes from hindbrain versus liver as measured by RNA-seq
#'
#' This dataset contains a dataframe of the top differentially expressed genes in the hindbrain versus liver as measured by RNA-seq (both genes that go up and those that go down). The gene names are the rownames, and the first column is the 'stat' column from DESeq2. Data was downloaded from ENCODE.
#'
#' \itemize{
#' \item gene_list_dataframe
#' }
#'
#' @docType data
#' @keywords datasets
#' @name gene_list_dataframe
#' @usage data(gene_list_dataframe)
#' @return A dataframe of top differentially expressed genes from hindbrain versus liver as measured by RNA-seq/
NULL

#' Character vector of the top differentially expressed genes from hindbrain versus liver as measured by RNA-seq
#'
#' This dataset contains a character vector of the top differenetially expressed genes in the hindbrain versus liver as measured by RNA-seq (both genes that go up and those that go down). Data was downloaded from ENCODE.
#'
#' \itemize{
#' \item gene_list_character
#' }
#'
#' @docType data
#' @keywords datasets
#' @name gene_list_character
#' @usage data(gene_list_character)
#' @return A character vector of the top differenetially expressed genes in the hindbrain versus liver as measured by RNA-seq/
NULL

#' GRangesList of the top 5000 H3K27ac peaks from hindbrain and liver downloaded from ENCODE
#'
#' This dataset contains a GRangesList of the H3K27ac peaks in either the hindbrain or the liver with the highest signal. Data was downloaded from ENCODE.
#'
#' \itemize{
#' \item K27ac_GRlist_hind_liver_top5000
#' }
#'
#' @docType data
#' @keywords datasets
#' @name K27ac_GRlist_hind_liver_top5000
#' @usage data(K27ac_GRlist_hind_liver_top5000)
#' @return A GRangesList of the top 5000 H3K27ac peaks from hindbrain and liver downloaded from ENCODE/
NULL


