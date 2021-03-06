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

  
  if (missing(i) && missing(j) && missing(k))
    return(x)
  

  if (!missing(i)) {
    if (is.character(i)) {
      fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
      i <- .profileplyr.charbound(i, rownames(x), fmt)
    }
    ii <- as.vector(i)
    ans_elementMetadata <- x@elementMetadata[i, , drop=FALSE]
    ans_rowRanges <- x@rowRanges[i]
    if(class(x@assays)[1] == "ShallowSimpleListAssays"){
        ans_assays <- x@assays[ii,]$data
    }
    else {
        ans_assays <- x@assays[ii,]@data
    }
    
    
    x  <- profileplyr_Dataset(ans_assays,ans_rowRanges,
                                         sampleData(x),sampleParams(x), params(x))

    
  }
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- .profileplyr.charbound(j, colnames(x), fmt)
    }
    ans_colData <- x@colData[j, , drop=FALSE]
    jj <- as.vector(j)
    if(class(x@assays)[1] == "ShallowSimpleListAssays"){
        ans_assays <- x@assays[,jj]$data
    }
    else {
        ans_assays <- x@assays[,jj]@data
    }
    
    x  <- profileplyr_Dataset(ans_assays,rowRanges(x),
                              sampleData(x),sampleParams(x),params(x))

    
  }
  if (!missing(k)) {
    if (is.character(k)) {
      fmt <- paste0("<", class(x), ">[,k] index out of bounds: %s")
      k <- .profileplyr.charbound(k, assayNames(x), fmt)
    }
    metaData <- x@metadata
    sampleDataTmp <- sampleData(x)[k, , drop=FALSE]
    sampleParamsTmp <- sampleParams(x)[k, , drop=FALSE]
    kk <- as.vector(k)
    if(class(x@assays)[1] == "ShallowSimpleListAssays"){
        ans_assays <- x@assays$data[kk]
    }
    else {
        ans_assays <- x@assays@data[kk]
    }

    x  <- profileplyr_Dataset(ans_assays,rowRanges(x),
                              sampleDataTmp,sampleParamsTmp,params(x))

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
                   if(is.null(rownames(sampleData(object)))) stop("Replacement sampleData must have rownames")

                   object@sampleData <- value
                   
                   names(assays(object)) <- rownames(value)
                   return(object)
                 })


sampleParams <- function (x)
          {
            return(x@sampleParams)
          }

#' Retrieve and set parameters in profileplyr object
#' @rdname params
#' @param object A profileplyr object
#' @return  A list containing parameters for profileplyr object.
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' params(object)
#' @export
params <- function (object)
{
  return(object@params)
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
                     params=params(x),
                     sampleParams=newSampleParams,
                     sampleData=newSampleData)
            metadata(x) <- metaX
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
#' @return A character vector of the top differentially expressed genes in the hindbrain versus liver as measured by RNA-seq/
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


