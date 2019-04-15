
#' Export and import profileplyr from/to deeptools
#'
#' Export and Import files
#'
#' @rdname export_deepToolsMat
#' @param object A profileplyr object
#' @param con Connection to write/read deeptools data to/from.
#' @param decreasing If object@params$mcolToOrderBy has been set and not NULL, then the ranges will be ordered by the column indicated in this slot of the metadata. By default, the order will be increasing for the factor or numeric value. For decreasing order, choose decreasing = TRUE.
#' @param overwrite Logical specifying whether to overwite output if it exists.
#' @importFrom R.utils gzip
#' @importFrom rjson toJSON
#' @importFrom dplyr mutate select vars
#' @importFrom tidyr gather separate
#' @importFrom magrittr "%>%"
#' @importClassesFrom methods  ANY character list missing
#' @importClassesFrom GenomicRanges GenomicRanges_OR_GRangesList
#' @importClassesFrom S4Vectors character_OR_NULL DataFrame
#' @importMethodsFrom GenomicFeatures as.list
#' @importMethodsFrom ComplexHeatmap draw
#' @importMethodsFrom GenomicRanges as.data.frame as.factor duplicated end findOverlaps match order score seqnames start strand
#' @importMethodsFrom IRanges as.matrix cbind colnames "colnames<-" diff gsub lapply ncol nrow paste quantile rbind rownames "rownames<-" subsetByOverlaps table unique which
#' @importMethodsFrom S4Vectors "%in%" mcols "mcols<-" do.call expand.grid grep grepl head levels metadata "metadata<-" Reduce subset t tail
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap anno_summary Heatmap HeatmapAnnotation
#' @importFrom dplyr mutate select
#' @importFrom EnrichedHeatmap anno_enriched EnrichedHeatmap
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom grid gpar unit
#' @importFrom IRanges IRanges
#' @importFrom magrittr "%>%"
#' @importFrom methods as is new
#' @importFrom pheatmap pheatmap
#' @importFrom R.utils gzip
#' @importFrom rGREAT submitGreatJob
#' @importFrom rjson fromJSON
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom soGGi regionPlot
#' @importFrom stats cutree

#' @importFrom methods as new
#' @importFrom utils combn read.delim write.table
#' @import SummarizedExperiment BiocGenerics
#' @export
#' @docType methods
#' @return  The path to deepTools matrix file
setGeneric("export_deepToolsMat", function(object="profileplyr",con="character",decreasing = "logical",overwrite = "logical") standardGeneric("export_deepToolsMat"))

#' @export
#' @describeIn export_deepToolsMat Export and import profileplyr from/to deeptools
setMethod("export_deepToolsMat", signature="profileplyr",
                                           function(object,con, decreasing = FALSE,overwrite=FALSE){
  con_prezip <- gsub("\\.gz","",con)
  con <- paste0(con_prezip,".gz")

  rowGroupsInUse <- params(object)$rowGroupsInUse
  mcolToOrderBy <- params(object)$mcolToOrderBy
  if(is.null(mcolToOrderBy)){
    order <- order(mcols(object)[, rowGroupsInUse])
  } else {
    
    order <- order(mcols(object)[, rowGroupsInUse],
                   mcols(object)[, mcolToOrderBy],
                   decreasing = decreasing)
  }
  
  object <- object[order, ]
  names(object) <- NULL
  
  groupInfo <- getGroupInfoFromObject(object)

  group_boundaries <- groupInfo$group_boundaries
  group_labels <- groupInfo$group_labels
  
  groupsL <- list(group_boundaries=group_boundaries,group_labels=group_labels)
  
  
  sbL <- list(sample_boundaries=assays(object) %>% 
                  lapply(ncol) %>% 
                  unlist %>% 
                  unname %>% 
                  cumsum %>% 
                  c(0,.))
  
  
  perSampleDPmetrics <- params(object)$perSampleDPParams
  
  perSampleDP <- sampleData(object) %>% 
                  as.data.frame %>% 
                  dplyr::select(!!perSampleDPmetrics) %>% 
                  as.list() %>% 
                  lapply(as.vector)
  
  perRun <- sampleData(object) %>% 
    as.data.frame %>% 
    dplyr::select(-!!perSampleDPmetrics) %>% lapply(unlist) %>% 
    lapply(unique) %>% lapply(as.vector)
  
  sc <- c(perSampleDP,perRun)
  names(sc) <- gsub("\\."," ",names(sc))
  
  if(nrow(sampleData(object)) == 1){
    toModNull <- sc[unlist(lapply(sc,is.null))]
    toMod <- sc[!unlist(lapply(sc,is.logical)) & lengths(sc) ==1]
    toModSC <- sc[unlist(lapply(sc,is.logical)) & lengths(sc) ==1]
    toGL <- sc[names(sc) =="group_labels"]
    sc <- c(lapply(toMod,function(x)unname(list(x))),toModSC,toModNull,toGL) 
  }
  sc <- c(sc,groupsL,sbL)
  
  if(file.exists(con) & !overwrite){
    stop(paste0(con_prezip," already exists. To overwrite file please set \"overwrite=TRUE\""))
  }
  
  sc %>% 
    rjson::toJSON(.) %>% 
      as.character %>%
      paste0("@",.) %>%
      write.table(con_prezip,row.names=FALSE,sep="",quote=FALSE,col.names=FALSE)



  scoreMat <- do.call(cbind,
                      as.list(SummarizedExperiment::assays(object)))


  as.data.frame(rowRanges(object)) %>% 
    mutate(score=".") %>% 
    select(seqnames,start,end,names,score,strand) %>%
    mutate(strand=gsub("\\*",".",strand)) %>%
    cbind(scoreMat) %>%
    write.table(con_prezip,row.names=FALSE,sep="\t",quote=FALSE,col.names=FALSE,append=TRUE)
  

  gzip(con_prezip,con,overwrite=TRUE,remove=FALSE)
  if(con_prezip !=con) unlink(con_prezip)
  
  return(con_prezip)
})

#' Import
#'
#'
#' @rdname export_deepToolsMat
#' @details A profileplyr object
#' @return A profileplyr object
#' @examples
#' 
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' export_deepToolsMat(object,file.path(tempdir(),"ATAC_Example.MAT"))
#' @importFrom rjson fromJSON
#' @export 
import_deepToolsMat <- function(con){
  myTempH <- scan(con,nlines = 1,what = "character",sep="\n")
  myTemp <- read.delim(con,sep="\t",comment.char="@",header=FALSE)
  myTempM <- as.matrix(myTemp[,-seq(6)])
  
  if(sum(is.na(myTempM)) > 0) {
    myTempM[is.na(myTempM)] <- 0
    message("Matrix contained 'NA' values, these will be replaced with zeros")
  }
  
  myTempGR <- GRanges(myTemp[,1],IRanges(myTemp[,2],myTemp[,3]),
                      names=myTemp[,4],
                      score=myTemp[,5],
                      strand=gsub("\\.","*",myTemp[,6]))
  
  info <- fromJSON(gsub("@","",myTempH))
  bounds <- sample_boundaries <- info$sample_boundaries
  sample_labels <- info$sample_labels
  
  sample_boundariesToFilter <- info$sample_boundaries
  sample_Starts <- sample_boundaries[-length(sample_boundaries)]
  sample_Ends <- sample_boundaries[-1]  
  sample_DataCol <- list()
  for(i in seq_along(sample_Starts)){
    sample_DataCol[[i]] <- seq(sample_Starts[i]+1,sample_Ends[i])
  }
  names(sample_DataCol) <- sample_labels
  
  myTempM_L <- lapply(sample_DataCol,function(x)myTempM[,x])
  myTempM_L <- lapply(myTempM_L,function(x){
    colnames(x) <- NULL
    x}
  )


  info_for_sampleData <- info[!(names(info) %in% c("group_labels", "group_boundaries"))]
  
  info_for_sampleData$`ref point` <- unlist(info_for_sampleData$`ref point`)
  if(length(sample_labels) > 1){

  sampleData <- DataFrame(as.data.frame(info_for_sampleData[lengths(info_for_sampleData) == length(sample_labels)]) %>%
                            cbind(info_for_sampleData %>% .[lengths(.) == 1] %>% lapply(rep,length(sample_labels)) %>% as.data.frame
                            ),
                          row.names = sample_labels)

  }else{
    sampleData <- DataFrame(as.data.frame(info_for_sampleData %>% .[lengths(.) == 1]
                              ),
                            row.names = sample_labels) 
  }
  if (is.null(unlist(info$`ref point`))){
    sampleData$ref.point <- replicate(nrow(sampleData),NULL)
  }
  sampleData$max.threshold <- sampleData$min.threshold <- replicate(nrow(sampleData),NULL)
  myTempGR$dpGroup <- factor(
    rep(info$group_labels,
        times=diff(info$group_boundaries) # Doug changed this too, before was 'each=diff(info$group_boundaries)' and this threw a warning that only first elsemnt was used. I think if you have different group sizes this wouldn't work, and 'times' fixes that
    ),
    levels = info$group_labels
  )
  
  perSampleDPParams <- info[lengths(info) == length(info$sample_labels)] %>% .[!(names(.) %in% c("group_labels", "group_boundaries"))] %>%
    names %>% make.names
  rowGroupsInUse <- "dpGroup"
  params <- list(perSampleDPParams=perSampleDPParams,rowGroupsInUse=rowGroupsInUse,mcolToOrderBy=rowGroupsInUse)
  
  proplyDataset <- profileplyr_Dataset(myTempM_L,myTempGR,sampleData,sampleParam=sampleData,params=params)
  return(proplyDataset)
}

#' Cluster Ranges
#'
#' Cluster the ranges in a deepTools object based on signal within each range
#' 
#' @docType methods
#' @name clusterRanges
#' @rdname clusterRanges
#' @param object A profileplyr object
#' @param fun The function used to summarize the ranges (e.g. rowMeans or rowMax)
#' @param scaleRows If TRUE, the rows matrix containing the signal in each bin that is used as the input for clustering will be scaled (as specified by pheatmap)
#' @param kmeans_k The number of kmeans groups used for clustering
#' @param clustering_callback Clustering callback function to be passed to pheatmap
#' @param clustering_distance_rows distance measure used in clustering rows. Possible values are "correlation" for Pearson correlation and all the distances supported by dist, such as "euclidean", etc. If the value is none of the above it is assumed that a distance matrix is provided.
#' @param cluster_method clustering method used. Accepts the same values as hclust
#' @param cutree_rows The number of clusters for hierarchical clustering
#' @param silent Whether or not a heatmap (from pheatmap) is shown with the output. This will not change what is returned with the function as it will always be a profileplyr object. If silent = FALSE, the heatmap will be shown which may be helpful in quick evaluation of varying numbers of clusters before proceeding with downstream analysis. The default is silent = TRUE, meaning no heatmap will be shown. 
#' @param show_rownames for any heatmaps printed while running this function, set to TRUE if rownames should be displayed. Default is FALSE. 
#' @details tbd
#' @return A profileplyr object
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' 
#' # k-means clustering
#' clusterRanges(object, fun = rowMeans, kmeans_k = 3)
#'
#' # hierarchical clustering, print heatmap, yet still return profileplyr object
#' clusterRanges(object, fun = rowMeans, cutree_rows = 3, silent = FALSE)
#'  
#' @export
#' @import IRanges
#' @importFrom pheatmap pheatmap
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @rdname clusterRanges
setGeneric("clusterRanges", function(object="profileplyr",fun="function",scaleRows="logical",
                                     kmeans_k="integer",clustering_callback="function",clustering_distance_rows="ANY",
                                     cluster_method="function",cutree_rows="integer",silent="logical",show_rownames="logical"
) standardGeneric("clusterRanges"))

#' @export
#' @describeIn clusterRanges Cluster Ranges
setMethod("clusterRanges", signature="profileplyr",
          function(object, fun = rowMeans, scaleRows = TRUE, kmeans_k = NULL, clustering_callback = function(x, ...){return(x)}, clustering_distance_rows = "euclidean", cluster_method = "complete", cutree_rows = NULL, silent = TRUE, show_rownames = FALSE) {

  
  if (length(assays(object)) < 2){
    stop("The profileplyr object must have at least two samples (length(assays(object)) > 1) for clustering")
  }
  # summarize adds the group name to the rownames, which I don't want to do as it messes up the separation required later on to build GRanges, but summarize is only a couple lines of code, so just include releveant ones here
  range_summ <- lapply(assays(object), fun)
  range_summ <- as.matrix(do.call(cbind,range_summ))

  colnames(range_summ) <- rownames(sampleData(object))
  
  if (scaleRows == TRUE) {
    scale <- "row"
  } else {
    scale <- "none"
  }
  
  if (is.null(kmeans_k) & is.null(cutree_rows)){
    res = pheatmap(range_summ, scale = scale, clustering_distance_rows = clustering_distance_rows, cluster_method = cluster_method, show_rownames = show_rownames, silent = FALSE, clustering_callback = clustering_callback)
    message("No 'kmeans_k' or 'cutree_rows' arguments specified. profileplyr object will be returned new column with hierarchical order from hclust")
    rowRanges(object)$hierarchical_order <- order(res$tree_row$order)
    return(object)
    }
  
  if (!(is.null(kmeans_k)) & !(is.null(cutree_rows))){
    message("Values for both 'kmeans_k' or 'cutree_rows' arguments were specified, kmeans will be used. To choose one method, change the value of the method not being used to NULL")
  }
  
  if (scaleRows == TRUE) {
    range_summ <- t(scale(t(range_summ)))
  }
  
  if (!(is.null(kmeans_k))) {
    message("K means clustering used.")
    res <- pheatmap(range_summ, scale = scale, kmeans_k = kmeans_k, silent = silent)
    cluster <- res$kmeans$cluster %>%
      data.frame(row.names = names(.), kmeans_cluster = .)
    rowRanges(object)$cluster <- cluster[,1]
  }else if (!(is.null(cutree_rows))) {
    if(identical(clustering_callback,(function(x, ...){return(x)}))){
      message("Hierarchical clustering used. It is advised to avoid this option with large matrices as the clustering can take a long time. Kmeans is more suitable for large matrices.")
    }
    if(!(identical(clustering_callback,(function(x, ...){return(x)})))){
      message("Hierarchical clustering performed using clustering method input in the 'callback_clustering' argument.  It is advised to avoid this option with large matrices as the clustering can take a long time. Kmeans is more suitable for large matrices.")
    }
    res <- pheatmap(range_summ, scale = scale, clustering_distance_rows = clustering_distance_rows, cluster_method = cluster_method, silent = silent, cutree_rows = cutree_rows, show_rownames = show_rownames, clustering_callback = clustering_callback)
    cluster <- cutree(res$tree_row, k = cutree_rows) %>%
      as.data.frame(.)
    rowRanges(object)$cluster <- cluster[,1]
    rowRanges(object)$hierarchical_order <- order(res$tree_row$order)
  }
  
  object@params$rowGroupsInUse <- "cluster"  # set the group to be filtered by as the overlap column, this will be used by export function to make the groups and also wil lbe labels for heatmap
  message("A column has been added to the range metadata with the column name 'cluster', and the 'rowGroupsInUse' has been set to this column.")         
  
  rowRanges(object)$cluster <- ordered(rowRanges(object)$cluster, 
                                                 levels = rowRanges(object)$cluster[!duplicated(rowRanges(object)$cluster)])
  
  colnames(mcols(object)) <- make.unique(colnames(mcols(object)))
  return(object)
}
)


#' Annotate profileplyr ranges to genes using rGREAT or ChIPseeker
#'
#' The ranges from the deepTools matrix will be subset based on whether they overlap with specified annotated regions related to a user defined gene list.
#' @docType methods
#' @name annotateRanges_great
#' @rdname annotateRanges
#' @param species GREAT accepts "hg19", "mm10", "mm9", "danRer7" (zebrafish)
#' @param ... pass to \code{\link[rGREAT]{submitGreatJob}} or \code{\link[ChIPseeker]{annotatePeak}}
#' @details tbd
#' @return A profileplyr object
#' @examples
#' 
#' library(SummarizedExperiment)
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example)
#' object <- object[1:2, , ] 
#' 
#' # annotate ranges with genes using ChIPseeker 
#' # (NOTE: can choose subset of annotations with 'annotation_subset' argument)
#' 
#' annotateRanges(object, TxDb = "mm10")
#' 
#' # annotate ranges with genes using GREAT with follwoing command (not run):
#' # annotateRanges_great(object, species = "mm10")
#' 
#' 
#' 
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Mmusculus.UCSC.mm9.knownGene TxDb.Mmusculus.UCSC.mm10.knownGene org.Hs.eg.db org.Mm.eg.db ChIPseeker rGREAT GenomicFeatures 
#' @rdname annotateRanges
setGeneric("annotateRanges_great", function(object="profileplyr",species="character",...)standardGeneric("annotateRanges_great"))


#' @describeIn annotateRanges Annotate profileplyr ranges to genes using rGREAT or ChIPseeker
#' @export
setMethod("annotateRanges_great", signature(object="profileplyr"),function(object, species, ...) {

  great <- submitGreatJob(rowRanges(object), request_interval = 0, species = species, ...)
  genomic_regions <- plotRegionGeneAssociationGraphs(great)
  colnames(mcols(genomic_regions)) <- c("SYMBOL", "distanceToTSS")
  
  hits <- findOverlaps(genomic_regions, rowRanges(object)) # get the indecies of both subject and query that overlap
  hits <- hits[!duplicated(queryHits(hits))]
  object <- object[subjectHits(hits)]
  mcols(object) <-  c(mcols(object), mcols(genomic_regions))
  
  colnames(mcols(object)) <- make.unique(colnames(mcols(object)))
  return(object)

})





#' Annotate profileplyr ranges to genes using rGREAT or ChIPseeker
#'
#' The ranges from the deepTools matrix will be subset based on whether they overlap with specified annotated regions (using ChIPseeker) throughout the genome
#' @docType methods
#' @name annotateRanges
#' @rdname annotateRanges
#' @param object A profileplyr object
#' @param annotation_subset If specific annotations (from ChIPseeker package) are desired, specify them here in a character vector. Can be one or any combination of "Promoter", "Exon", "Intron", "Downstream", "Distal Intergenic", "3p UTR", or "5p UTR". This argument is optional and all annotation types will be included if argument is left out.
#' @param TxDb TxDb object, a character sting that is a path to a GTF file, or character string indicating genome if one of the following - "hg19", "hg38", "mm9", "mm10".
#' @param tssRegion This needs to be a vector of two numbers that will define promoter regions. The first number must be negative, while the second number must be positive. Default values are  c(-3000, 3000) - SHOULD WE CHANGE THIS, SEEMS BIG!)
#' @param changeGroupToAnnotation If the grouping should be changed to the annotations (typically when the ranges will be exported for visualization based on this annotation), this should be TRUE. The default if FALSE, which will keep the grouping that existed before annotating the object. This is typical if the output will be used for finding overlaps with gene lists in the 'groupBy' function.
#' @param heatmap_grouping Only relevant if 'keepAnnotationAsGroup' is set to TRUE. This argument needs to be either "group", or "annotation". This will determine how the ranges are grouped in the resulting object. Default is heatmap_grouping = "Group". If there are no groups in the deepTools matrix that was used in the function, this argument is unnecessary
#' @details tbd
#' @return A profileplyr object
#' @rdname annotateRanges 
#' @export
setGeneric("annotateRanges", function(object="profileplyr",annotation_subset = "character", TxDb, tssRegion = "numeric", changeGroupToAnnotation = "logical", heatmap_grouping = "character", ...)standardGeneric("annotateRanges"))

#' @describeIn annotateRanges Annotate profileplyr ranges to genes using rGREAT or ChIPseeker
#' @export
setMethod("annotateRanges", signature(object="profileplyr"),function(object, annotation_subset = NULL, TxDb, tssRegion = c(-3000, 3000), changeGroupToAnnotation = FALSE, heatmap_grouping = "group", ...) {
  
  
  # which TxDb?
  if (is(TxDb,"TxDb")) {
    TxDb_object <- TxDb
  } else if (is(TxDb,"character") & file.exists(TxDb)) {
    TxDb_object <- makeTxDbFromGFF(TxDb)
  } else if(is(TxDb,"character") & !(file.exists(TxDb))) {
    if (TxDb == "hg19") {
      TxDb_object <- TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (TxDb == "hg38") {
      TxDb_object <- TxDb.Hsapiens.UCSC.hg38.knownGene
    } else if (TxDb == "mm9") {
      TxDb_object <- TxDb.Mmusculus.UCSC.mm9.knownGene
    } else if (TxDb == "mm10") {
      TxDb_object <- TxDb.Mmusculus.UCSC.mm10.knownGene
    }
  }
  
  
  # which species for annotation?
  if (TxDb == "hg19" | TxDb == "hg38") {
    orgAnn <- "org.Hs.eg.db"
  } else if (TxDb == "mm9" | TxDb == "mm10") {
    orgAnn <- "org.Mm.eg.db"
  } else {
    orgAnn <- NULL
  }
  
  peakanno <- annotatePeak(rowRanges(object), TxDb = TxDb_object, annoDb = orgAnn, tssRegion = tssRegion, ...)
  annotatedPeaksGR <- as.GRanges(peakanno)
  
  annotatedPeaksGR$range_combined <- paste(seqnames(annotatedPeaksGR), start(annotatedPeaksGR), end(annotatedPeaksGR), sep = "_")
  rowRanges(object)$range_combined <- paste(seqnames(rowRanges(object)), start(rowRanges(object)), end(rowRanges(object)), sep = "_")
  
  object <- object[rowRanges(object)$range_combined %in% annotatedPeaksGR$range_combined]
  mcols(object) <- mcols(annotatedPeaksGR)
  
  mcols(object) <-  subset(mcols(object), select=-c(range_combined))
  
  rowRanges(object)$annotation_short <- gsub("\\s*\\([^\\)]+\\)","",as.character(rowRanges(object)$annotation))
  
  # output matrix throws error in plotHeatmap call with the single quote, as it interprets as a real quote, so change to "p"
  rowRanges(object)$annotation_short <- gsub("3' UTR", "3p UTR",as.character(rowRanges(object)$annotation_short))
  rowRanges(object)$annotation_short <- gsub("5' UTR", "5p UTR",as.character(rowRanges(object)$annotation_short))
  annotation_order <- c("Promoter", "Exon", "Intron", "Downstream", "Distal Intergenic", "3p UTR", "5p UTR")
  rowRanges(object)$annotation_short <- ordered(rowRanges(object)$annotation_short, levels = annotation_order)

  # if the user wants a subset of groups, this will subset and also order the groups based on the order of input
  if (!(is.null(annotation_subset))){
    object <- object[rowRanges(object)$annotation_short %in% annotation_subset]
    rowRanges(object)$annotation_short <- ordered(rowRanges(object)$annotation_short, levels = annotation_subset)
  }
  
  if (changeGroupToAnnotation == TRUE){
    
    group <- params(object)$rowGroupsInUse
    group_label <- mcols(object)[colnames(mcols(object)) %in% group][,1]
    
    if(heatmap_grouping == "group") {
      object <- object[order(group_label, rowRanges(object)$annotation_short)] # order first based on group, then by annotation
      group_label_order <- mcols(object)[colnames(mcols(object)) %in% group][,1] # get the new group_label order based on the above ordering
      rowRanges(object)$annotate_group <- paste(group_label_order, rowRanges(object)$annotation_short, sep = "_") # make new column with combined grouping/annotation
    }
    
    if(heatmap_grouping == "annotation") {
      object <- object[order(rowRanges(object)$annotation_short, group_label)]
      group_label_order <- mcols(object)[colnames(mcols(object)) %in% group][,1]
      rowRanges(object)$annotate_group <- paste(rowRanges(object)$annotation_short, group_label_order, sep = "_")
    }
    object@params$rowGroupsInUse <- "annotate_group"
    message("A column has been added to the range metadata with the column name 'annotate_group' that combines the inherited groups with the annotations determined here. The 'rowGroupsInUse' has been set to this column.")      
  }
  
  colnames(mcols(object)) <- make.unique(colnames(mcols(object)))
  return(object)
  
})





#' Redundant code wrapped into subsetbyRangeOverlap() or subsetbyGeneListOverlap() functions
#'
#' 
#'
#' @rdname subset_GR_GL_common_top
#' @param object A profileplyr object
#' @param overlap hits object from subsetByOverlap function
#' @param input_names names of either the gene list of the granges that go into function 
#' @param type Either "GR" for subsetbyRangeOverlap() function or "GL" for subsetbyGeneListOverlap() function
#' @param separateDuplicated A logical argument, if TRUE then regions that overlap multiple inputs to 'GRanges' argument will be separated and made into their own group. All possible combinations of region overlaps will be tested, so it is not recommended to have more than 3 groups if this option is TRUE. If FALSE, then regions that overlap each individual 'GRanges' input will be in the output, and if one region overlaps multiple 'GRanges' inputs, then it will be duplicated in the output and will show up in the section for each group.
#' @details tbd
#' @return A list of profileplyr objects
#' 
#'  
subset_GR_GL_common_top <- function(object, overlap, input_names, type, separateDuplicated) {
  
  # this code will give us a logical for the rows that overlap nothing
  overlap_combined <- do.call(cbind, overlap)
  no_overlap <- rowSums(overlap_combined) == 0
  
  # make a matrix with a 0 or 1 saying whether there was an overlap with each set of ranges
  overlap_matrix <- matrix(nrow = length(overlap[[1]]), ncol = length(overlap))
  for(i in seq_along(overlap)){
    overlap_matrix[, i] <- as.numeric(overlap[[i]])
  }
  
  colnames(overlap_matrix) <- input_names
  
  rowData(object)$overlap_matrix <- overlap_matrix # add column in rowData containing this matrix
  
  # following code serves purpose of adding a column to rowData that will show the grouping, including all combinations of overlaps. So each range only has one assignment
  overlap_names <- vector()
  for (i in seq_len(nrow(overlap_matrix))){
    if (sum(overlap_matrix[i,]) == 0) {
      overlap_names[i] <- "no_overlap"
    } else {
      overlap_names[i] <- names(overlap_matrix[1,]) %>%
        .[as.logical(overlap_matrix[i,])] %>% #this will index the names the input GRanges based on each row of the matrix
        paste(., collapse = "_")
    }
  }
  rowRanges(object)$overlap_names <- overlap_names
  rowRanges(object)$overlap_names <- ordered(rowRanges(object)$overlap_names, 
                                             levels = rowRanges(object)$overlap_names[!duplicated(rowRanges(object)$overlap_names)])
  
  object_overlap <-  object[which(!no_overlap)]
  object_no_overlap <- object[no_overlap]
  
  colnames(mcols(object_overlap))[colnames(mcols(object_overlap)) %in% "overlap_names"] <- paste0(type, "_overlap_names")
  colnames(mcols(object_no_overlap))[colnames(mcols(object_no_overlap)) %in% "overlap_names"] <- paste0(type, "_overlap_names")
  
  if (separateDuplicated == TRUE) {
    object_overlap@params$rowGroupsInUse <- paste0(type, "_overlap_names")
    object_no_overlap@params$rowGroupsInUse <- paste0(type, "_overlap_names")
    message(paste0("A column has been added to the range metadata with the column name '", type, "_overlap_names' that specifies the GRanges each range overlaps with, but the inherited groups are not included."))      
  }
  
  if (separateDuplicated == FALSE) {
    
    overlap_matrix_df <- as.data.frame(overlap_matrix)
    overlap_list <- lapply(overlap_matrix_df, function(x) object[which(as.logical(x))])
    
    # make a vector with the names of each GRanges with the length of each overlap profileplyr object
    overlap_nosep_names <- vector()
    for(i in seq_along(overlap_list)){
      temp <- rep(colnames(overlap_matrix)[i], length(overlap_list[[i]]))
      overlap_nosep_names <- c(overlap_nosep_names, temp)
    }
    
    object_overlap <- do.call(rbind, overlap_list)
    mcols(object_overlap)$overlap_nosep_names <- overlap_nosep_names
    rowRanges(object_overlap)$overlap_nosep_names <- ordered(rowRanges(object_overlap)$overlap_nosep_names, 
                                                             levels = rowRanges(object_overlap)$overlap_nosep_names[!duplicated(rowRanges(object_overlap)$overlap_nosep_names)])
    
    object_overlap@params$rowGroupsInUse <- paste0(type, "_overlap_nosep_names")
    object_no_overlap@params$rowGroupsInUse <- paste0(type, "_overlap_nosep_names")
    message(paste0("A column has been added to the range metadata with the column name '", type, "_overlap_nosep_names' that specifies the GRanges each range overlaps with, but the inherited groups are not included."))      
    
    colnames(mcols(object_overlap))[colnames(mcols(object_overlap)) %in% "overlap_names"] <- paste0(type, "_overlap_names")
    colnames(mcols(object_overlap))[colnames(mcols(object_overlap)) %in% "overlap_nosep_names"] <- paste0(type, "_overlap_nosep_names") 
    
    mcols(object_no_overlap)$overlap_nosep_names <- "no_overlap"
    colnames(mcols(object_no_overlap))[colnames(mcols(object_no_overlap)) %in% "overlap_nosep_names"] <- paste0(type, "_overlap_nosep_names")
    
  }
  return(list(object_overlap, object_no_overlap))
}

#' Redundant code for inheriting grouping wrapped into subsetbyRangeOverlap() or subsetbyGeneListOverlap() functions
#'
#' 
#'
#' @rdname inherit_group_function
#' @param object A profileplyr object
#' @param rowGroupsInUse_input the inherited rowGroupsInUse
#' @param type Either "GR" for subsetbyRangeOverlap() function or "GL" for subsetbyGeneListOverlap() function
#' @param separateDuplicated A logical argument, if TRUE then regions that overlap multiple inputs to 'GRanges' argument will be separated and made into their own group. All possible combinations of region overlaps will be tested, so it is not recommended to have more than 3 groups if this option is TRUE. If FALSE, then regions that overlap each individual 'GRanges' input will be in the output, and if one region overlaps multiple 'GRanges' inputs, then it will be duplicated in the output and will show up in the section for each group.
#' @details tbd
#' @return A profileplyr object
#' 
inherit_group_function <- function(object, rowGroupsInUse_input, type, separateDuplicated) {
  
  inherit_groups <- mcols(object)[colnames(mcols(object)) %in% rowGroupsInUse_input] %>%
    .[,1] %>%
    as.character()

# create a column in the non overlapping GRanges that has the group and non-overlapping indication
if (separateDuplicated == TRUE){
  rowRanges(object)$group_and_overlap <- paste(inherit_groups,
                                                  as.character(mcols(object)[,colnames(mcols(object)) %in% paste0(type, "_overlap_names")]),
                                                  sep = "_")
  colnames(mcols(object))[colnames(mcols(object)) %in% "group_and_overlap"] <- paste0(type, "_group_and_overlap")
  object@params$rowGroupsInUse <- paste0(type, "_group_and_overlap")
  message(paste0("A column has been added to the range metadata with the column name '", type, "_group_and_overlap' that specifies the GRanges each range overlaps with, but the inherited groups are not included."))      
  
}
  
if (separateDuplicated == FALSE){
  rowRanges(object)$group_and_overlap_nosep <- paste(inherit_groups,
                                               as.character(mcols(object)[,colnames(mcols(object)) %in% paste0(type, "_overlap_nosep_names")]),
                                               sep = "_")
  colnames(mcols(object))[colnames(mcols(object)) %in% "group_and_overlap_nosep"] <- paste0(type, "_group_and_overlap_nosep") 
  object@params$rowGroupsInUse <- paste0(type, "_group_and_overlap_nosep")
  message(paste0("A column has been added to the range metadata with the column name '", type, "_group_and_overlap_nosep' that specifies the GRanges each range overlaps with, but the inherited groups are not included."))      
}
return(object)
}

#' Subset ranges based on overlap with a GRanges object
#'
#' The ranges from the deepTools matrix will be subset based on whether they overlap with user defined ranges
#'
#' @rdname subsetbyRangeOverlap
#' @param object A profileplyr object
#' @param group The regions by which to subset the deepTools matrix. This must be either a single GRanges object, or a GRangesList. Combinations of bed file paths and GRanges objects are not accepted, Import BED files as GRanges with rtracklayer import.bed() function.
#' @param GRanges_names The names of the GRanges that were used for the "GRanges" argument. This will be used to label these groups in the construction of the resulting profileplyr object.
#' @param include_nonoverlapping A logical argument, if FALSE the regions from the original deepTools matrix that do not overlap with the  user defined regions will be left out of the returned profileplyr object.
#' @param separateDuplicated A logical argument, if TRUE then regions that overlap multiple inputs to 'GRanges' argument will be separated and made into their own group. All possible combinations of region overlaps will be tested, so it is not recommended to have more than 3 groups if this option is TRUE. If FALSE, then regions that overlap each individual 'GRanges' input will be in the output, and if one region overlaps multiple 'GRanges' inputs, then it will be duplicated in the output and will show up in the section for each group.
#' @param inherit_groups A logical whether that groups the exist in the profileplyr object in the 'object' argument should be included in the default grouping scheme for the output object of this function. The default is TRUE. If false, only the GRanges overlap annotation will be used for heatmap grouping.
#' @details tbd
#' @return A profileplyr object
#' @examples
#' # see the groupby function within profileplyr for examples
#'  
subsetbyRangeOverlap <- function(object, group, GRanges_names = NULL, include_nonoverlapping = FALSE, separateDuplicated = TRUE, inherit_groups = FALSE) {
  
  rowGroupsInUse_input <- params(object)$rowGroupsInUse
  region_list_GRanges <- GRangesList(group)
  
  if(!(is.null(names(region_list_GRanges)))){
    GRanges_names <- names(region_list_GRanges)
  }else if(is.null(names(region_list_GRanges)) & is.null(GRanges_names)){
    message("The argument 'GRanges_names' was not set so the GRanges will be given generic names. To name GRanges, set them using the 'GRanges_names' argument")
    GRanges_names <- vector()
    for(i in seq_along(region_list_GRanges)) {
      temp <- paste0("RegionSet_", i)
      GRanges_names <- c(GRanges_names, temp)
      names(region_list_GRanges) <- GRanges_names
    }
  }else {
    names(region_list_GRanges) <- GRanges_names
  }
  
  if (length(region_list_GRanges) != length(GRanges_names)) {
    stop("The number of region sets for overlap analysis does not match the number of names for those sets")
  }
  
  # this will give us a list of logical vectors that tell us if a range at that index overlaps with each GRanges
  overlap <- lapply(region_list_GRanges, function(x) rowRanges(object) %in% subsetByOverlaps(rowRanges(object), x)) # get the ranges that are in the input range files and overlap with the object ranges (from the deepTools matrix)
  
  input_names = GRanges_names
  type = "GR"
  
  object_list <- subset_GR_GL_common_top(object = object, 
                                         overlap = overlap, 
                                         input_names = input_names, 
                                         type = type, 
                                         separateDuplicated)
   
  object_overlap <- object_list[[1]]
  object_no_overlap <- object_list[[2]]
  
  # add metadata columns from the overlapping GRanges to the profile plyr object
  region_list_GRanges_unlist <- unlist(region_list_GRanges)
  
  hits <- findOverlaps(object_overlap, region_list_GRanges_unlist) # get the indecies of both subject and query that overlap
  uniqueInQuery <- subjectHits(hits) %>%
    .[!duplicated(queryHits(hits))] # this will give us indecies of the input GRanges for subsetting (query) that overlap the deepTools SE (subject).
  query_mcols <- mcols(region_list_GRanges_unlist)[uniqueInQuery, ] # get the metadata from the GRanges for subsetting (query) that corresponds to the ranges that overlap with the deepTools SE (subject). If multiple input regions overlap with the SE, metadata from only the first one will be caught here
  mcols(object_overlap) <-  c( mcols(object_overlap), query_mcols) # combine the metadata from the subsetting GRanges with the metadata in the SE DT
  
  if (include_nonoverlapping == FALSE) {
    object <- object_overlap
  }
  
  if (include_nonoverlapping == TRUE) {
    
    # get a vector of the breakdown of groups within the non-overlapping ranges to make the same 'group_overlap' column we have the the overlap ranges above
    
    
    
    
    # here I just populate the metadata columns in the non-overlapping regions with NAs
    no_overlap_mcols <- list()
   
     for (i in seq_along(query_mcols)) {
      no_overlap_mcols[[i]] <- rep(NA, length(rowRanges(object_no_overlap)))
    }
    no_overlap_mcols <- do.call("cbind",no_overlap_mcols) %>%  as.data.frame()
    colnames(no_overlap_mcols) <- colnames(query_mcols)
    mcols(object_no_overlap) <-  c( mcols(object_no_overlap), no_overlap_mcols)
    
    object <- rbind(object_overlap, object_no_overlap)
  }
  if (inherit_groups == TRUE) {
    
    object <- inherit_group_function(object, 
                                     rowGroupsInUse_input, 
                                     type,
                                     separateDuplicated)
  }
  
  colnames(mcols(object)) <- make.unique(colnames(mcols(object)))
  return(object)
}

#' Subset ranges based on overlap with lists of Gene sets
#'
#' The ranges from the deepTools matrix will be subset based on whether they overlap with user defined gene sets
#'
#' @rdname subsetbyGeneListOverlap
#' @param object A profileplyr object
#' @param group The regions by which to subset the deepTools matrix. This must be either a single GRanges object, or a GRangesList. Combinations of bed file paths and GRanges objects are not accepted, Import BED files as GRanges with rtracklayer import.bed() function.
#' @param include_nonoverlapping A logical argument, if FALSE the regions from the original deepTools matrix that do not overlap with the  user defined regions will be left out of the returned profileplyr object.
#' @param separateDuplicated A logical argument, if TRUE then regions that overlap multiple inputs to 'GRanges' argument will be separated and made into their own group. All possible combinations of region overlaps will be tested, so it is not recommended to have more than 3 groups if this option is TRUE. If FALSE, then regions that overlap each individual 'GRanges' input will be in the output, and if one region overlaps multiple 'GRanges' inputs, then it will be duplicated in the output and will show up in the section for each group.
#' @param inherit_groups A logical whether that groups the exist in the profileplyr object in the 'object' argument should be included in the default grouping scheme for the output object of this function. The default is TRUE. If false, only the gene list overlap annotation will be used for heatmap grouping.
#' @details tbd
#' @return A profileplyr object
#' @examples
#' # see the groupby function within profileplyr for examples
#' 

subsetbyGeneListOverlap <- function(object, group, include_nonoverlapping = FALSE, separateDuplicated = TRUE, inherit_groups = FALSE) {
  
  rowGroupsInUse_input <- params(object)$rowGroupsInUse
  
  if(is.null(names(group))){
    message("Input gene lists do not have names, they will be given generic names. To name gene list, set them before using this function with names(gene list) function")
    gene_list_names <- vector()
    for(i in seq_along(group)) {
      temp <- paste0("GeneSet_", i)
      gene_list_names <- c(gene_list_names, temp)
    }
  } else {
    gene_list_names <- names(group)
  }

  # check to see whether all are data frames
  # if not , then we will just use gene symbols and no column data
  df_check_list <- lapply(group, is.data.frame) %>% 
                                            unlist() %>%
                                              all()
  
  #this will first check is all elements are data frames and length of 1, meaning it will take extra columns becuase it can't be differnet form toher data frames as there are none
  # it will then check for the all data frame lists larger than 1 if all the column names are the same
  # if either of the above are TRUE this is recorded and will be used later on
  if(df_check_list & length(group) == 1){
    keep_columns <- TRUE
  }else if(df_check_list & length(group) > 1){
    colnames_list <- lapply(group, colnames)
    identical_colnames <- lapply(colnames_list[2:length(colnames_list)], FUN = identical, colnames_list[[1]]) %>%
                            unlist() %>%
                              all()
    if (identical_colnames) {
      keep_columns <- TRUE
    }
  }else{
    keep_columns <- FALSE
  }

    # only the gene symbols will be extracted for overlaps
    gene_list_vectorList <- vector(mode = "list", length = length(group))
    for (i in seq_along(group)){
      if(is.data.frame(group[[i]])){
        gene_list_vectorList[[i]] <- row.names(group[[i]])
      }else if (is.character(group[[i]])){
        gene_list_vectorList[[i]] <- group[[i]]
      }else{
        stop("elements of 'group' must be a character vector, or a data frame with the gene symbols as the rownames")
      }
    }
  
    overlap <- lapply(gene_list_vectorList, function(x) rowRanges(object)$SYMBOL %in% x)
    
    input_names = gene_list_names
    type = "GL"
    
    object_list <- subset_GR_GL_common_top(object = object, 
                                           overlap = overlap, 
                                           input_names = input_names, 
                                           type = type, 
                                           separateDuplicated)  
    
    object_overlap <- object_list[[1]]
    object_no_overlap <- object_list[[2]]
    
  if(keep_columns){
    temp <- group
    names(temp) <- NULL
    temp <- do.call(rbind, temp)
    
    index <- match(rowRanges(object_overlap)$SYMBOL, row.names(temp))
    overlap_in_gene <- temp[index, , drop = FALSE]
    mcols(object_overlap) <- cbind(mcols(object_overlap), overlap_in_gene)
  }
  
  if (include_nonoverlapping == FALSE) {
    object <- object_overlap
  }
  
  
  if (include_nonoverlapping == TRUE) {
    
    object <- rbind(object_overlap, object_no_overlap)
    
  }
    
    if (inherit_groups == TRUE) {
      
      object <- inherit_group_function(object, 
                                       rowGroupsInUse_input, 
                                       type,
                                       separateDuplicated)
    }
  colnames(mcols(object)) <- make.unique(colnames(mcols(object)))
  return(object)
  
}

#' summarize the rows of a deepTools matrix
#'
#' summarize the rows of a deepTools matrix
#' @docType methods
#' @name summarize
#' @rdname summarize
#' @param object A profileplyr object
#' @param fun the function used to summarize the ranges (e.g. rowMeans or rowMax)
#' @param output Must be either "matrix", "long", or "object".
#' @param keep_all_mcols if output is 'long' and this is set to TRUE, then all metadata columns in the rowRanges will be included in the output. If FALSE (default value), then only the column indicated in the 'rowGroupsInUse' slot of the metadata will be included in the output dataframe. 
#' @details Takes a SE object and outputs a summarized experiment object with a matrix containing ranges as rows and each sample having one column with summary statistic
#' @return If output="matrix" returns a matrix, if output="long" returns a data.frame in long format,  if output="long" returns a SummarizedExperiment object
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example)
#' 
#' # output matrix (can be used to make a heatmap)
#' 
#' object_sumMat <- summarize(object, fun = rowMeans, output = "matrix") 
#' 
#' # output long dataframe for ggplot
#' 
#' object_long <- summarize(object, fun = rowMeans, output = "long") 
#' object_long[1:3, ]
#' 
#' library(ggplot2)
#' ggplot(object_long, aes(x = Sample, y = log(Signal))) + geom_boxplot()
#' 
#' # output profileplyr object containing summarized matrix
#' 
#' summarize(object, fun = rowMeans, output = "object")
#' @export 
setGeneric("summarize", function(object="profileplyr",fun = "function", output = "character", keep_all_mcols = "logical")standardGeneric("summarize"))

#' @describeIn  summarize summarize the rows of a deepTools matrix
#' @export
setMethod("summarize", signature(object="profileplyr"), function(object, fun, output, keep_all_mcols = FALSE){
  summ <- lapply(assays(object), fun)
  summ_mat <- as.matrix(do.call(cbind,summ))
  rowGroupsInUse <- params(object)$rowGroupsInUse
  GroupNames <- rowData(object)[,colnames(rowData(object)) %in% rowGroupsInUse]
  
  colnames(summ_mat) <- rownames(sampleData(object))
  
  if (output == "matrix") {
    rownames(summ_mat) <- paste(seqnames(rowRanges(object)), start(rowRanges(object)), end(rowRanges(object)), GroupNames, sep = "_") # not sure how critrical this is, but included group name in case ranges are in multiple groups
    return(summ_mat)
  }
  
  ####  this can be an output of the plotRanges function, or almost should be a function on its own. I feel like it might get buried here as it doesn't quite fit with the matrix/SE output
  if (output == "long") {
    summ_long <- as.data.frame(summ_mat)
    if (keep_all_mcols == TRUE){
      summ_long <- cbind(summ_long, mcols(rowRanges(object))) %>%
                                      data.frame(check.names = FALSE)
    } else {
      summ_long[rowGroupsInUse] <- GroupNames
    }
    summ_long$combined_ranges <- paste(seqnames(rowRanges(object)), start(rowRanges(object)), end(rowRanges(object)), sep = "_") # NOTE: don't include group names here as we have that in a separate column
    summ_long <- tidyr::gather(summ_long, key = "Sample", value = "Signal", seq_len(ncol(summ_mat)))
    summ_long$Sample <- ordered(summ_long$Sample, levels = rownames(sampleData(object)))
    return(summ_long)
  }
  
  if (output == "object") {
    
    objectToReturn <- SummarizedExperiment(summ_mat,rowRanges(object),
                                           colData=sampleData(object))
    metadata(objectToReturn) <- metadata(object)
    return(objectToReturn)
  }
})






#' group the rows and ranges of the profileplyr object
#'
#' group the rows and ranges of the profileplyr object
#' @docType methods
#' @name groupBy
#' @rdname groupBy
#' @param object  A profileplyr object
#' @param group How the ranges will be grouped. If this is a character string, then it must match a column name of the range metadata, and this column will be used for grouping of any exported deepTools matrix. If this is a GRanges, or GRangesList, then the ranges will be subset based on overlap with these GRanges. If this is a list, each element should contain a character vector of genes, and ranges will be subset based on overlap with these genes, as determined by the annotations made by annotateRanges() or annotateRanges_great() functions.
#' @param levels This will set the levels of the grouping column set by 'rowGroupsInUse' (if the grouping column is not a factor, it will be converted to one). If levels are not provided, they will remain unchanged if the grouping column was already a factor, or will use default leveling (e.g. alphabetical) if grouping column is not already a factor variable. 
#' @param GRanges_names Relevant for 'GRanges' mode. These are the names that will be assigned to the ranges that overlap each GRanges object
#' @param include_nonoverlapping Relevant for 'GRanges' mode. This should be indicated (default is TRUE). A logical argument, if FALSE the regions from the original deepTools matrix that do not overlap with the  user defined regions will be left out of the returned profileplyr object.
#' @param separateDuplicated Relevant for 'GRanges' mode. This should be indicated (default is FALSE). A logical argument, if TRUE then regions that overlap multiple inputs toe 'regions' argument will be separated and made into their own group. All possible combinations of region overlaps will be tested, so it is not recommended to have more than 3 groups if this option is TRUE. If FALSE, then regions that overlap each individual 'region' input will be in the output, and if one region overlaps multiple 'region' inputs, then it will be duplicated in the output and will show up in the section for each group.
#' @param inherit_groups A logical whether that groups the exist in the profileplyr object in the 'object' argument should be included in the default grouping scheme for the output object of this function. The default is TRUE. If false, only the GRanges or gene list overlap annotation will be used for heatmap grouping.
#' @details Takes a SE object and groups rows
#' @return  A profileplyr object with a summarized matrix, a matrix, or a long dataframe.
#' @examples
#' 
#' # group by gene list or list of data frames with genes as rownames
#' ## not shown here but see vignette for grouping by gene lists
#' 
#' # group by GRanges
#' 
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example)
#' data("K27ac_GRlist_hind_liver_top5000") # load pre-made GRanges
#' K27ac_groupByGR <- groupBy(object, group = K27ac_GRlist_hind_liver_top5000)
#' 
#' # switch rowGroupsInUse
#' 
#' switchGroup <- groupBy(K27ac_groupByGR, group = "GR_overlap_names")
#' metadata(switchGroup)$rowGroupsInUse
#' @export
setGeneric("groupBy", function(object="profileplyr",group="ANY", GRanges_names = "character", levels = "ANY", 
                               include_nonoverlapping = "logical", separateDuplicated = "logical", inherit_groups = "logical")standardGeneric("groupBy"))

#' @describeIn groupBy group the rows and ranges of the profileplyr object
#' @export
setMethod("groupBy", signature(object="profileplyr"),function(object, group, GRanges_names = NULL, levels = NULL, include_nonoverlapping = FALSE, separateDuplicated = TRUE, inherit_groups = FALSE){
  
  if(missing(group)){
    stop("Enter 'group' argument")
  }
  
  if(is(group,"character")){
    if(!(group %in% colnames(mcols(object)))) {
      stop("The 'group' argument is a character string, but does not match the name of any range metadata columns")
    }
    object@params$rowGroupsInUse <- group
    mcols(object)[, group] <- as.factor(mcols(object)[, group]) #make sure this is a factor so we can order it
    if (is.null(levels)){
      levels <- levels(rowData(object)[, group]) #if user doesnt set levels, we will just use the default ones (alphabetical)
    }
    mcols(object)[, group] <- ordered(mcols(object)[, group], levels = levels)

    return(object)
  }
  
  
  if(is(group,"GRanges") | is(group,"CompressedGRangesList")){
    return(subsetbyRangeOverlap(object, group = group, GRanges_names, include_nonoverlapping, separateDuplicated, inherit_groups))
  }
  
  if(is(group,"list")) {
    return(subsetbyGeneListOverlap(object, group = group, include_nonoverlapping, separateDuplicated, inherit_groups))
  }
 
})




#' choose the column by which to order the ranges by within each group
#'
#' choose the column by which to order the ranges by within each group
#' @docType methods
#' @name orderBy
#' @rdname orderBy
#' @param object  A profileplyr object
#' @param column Which column of mcols(proplyrObject) should be used for ordering the ranges
#' @details Takes a profileplyr object and orders the rows based on a user defined metadata column of rowRanges
#' @return  A profileplyr object
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' 
#' library(SummarizedExperiment)
#' cluster <- clusterRanges(object, fun = rowMeans, cutree_rows = 3)
#' cluster_order <- orderBy(cluster, column = "hierarchical_order")
#' metadata(cluster_order)$mcolToOrderBy
#'
#' @export
setGeneric("orderBy", function(object="profileplyr",column = "character")standardGeneric("orderBy"))

#' @describeIn orderBy choose the column by which to order the ranges by within each group
#' @export
setMethod("orderBy", signature(object="profileplyr"),function(object, column){

    if(!(column %in% colnames(mcols(object)))) {
      stop("The 'column' argument does not match the name of any range metadata columns")
    }
    # object@params$mcolToOrderBy <- column
    object@params$mcolToOrderBy <- column
    
    return(object)
})





#' export a profileplyr object to a list of matrices that can be used as an input for EnrichedHeatmap
#'
#' export a profileplyr object to a list of matrices that can be used as an input for EnrichedHeatmap
#'
#' @rdname convertToEnrichedHeatmapMat
#' @param object  A profileplyr object
#' @param sample_names A character vector that will set the names of the heatmap components that are generated from the profileplyr assays() matrices. This argument is optional, by default the names will be the name of the samples in the profileplyr object metadata(proplyrObject)$sampleData$sample_labels.
#' @details Takes a profileplyr object and converts all of the matrices in the assays() section of the object to matrices that can be used as an input for EnrichedHeatmap
#' @return A list of normalized matrices that can be used for generating visualizations with EnrichedHeatmap
#' @examples
#' 
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example) 
#' 
#' library(EnrichedHeatmap)
#' EH_mat <- convertToEnrichedHeatmapMat(object)
#' EnrichedHeatmap(EH_mat[[1]], name = names(EH_mat[1]), column_title = names(EH_mat[1]))
#' @import EnrichedHeatmap ComplexHeatmap grid circlize
#' @export
setGeneric("convertToEnrichedHeatmapMat", function(object="profileplyr",sample_names="character") standardGeneric("convertToEnrichedHeatmapMat"))

#' @describeIn convertToEnrichedHeatmapMat export a profileplyr object to a list of matrices that can be used as an input for EnrichedHeatmap
#' @export
setMethod("convertToEnrichedHeatmapMat", signature(object="profileplyr"),function(object, sample_names = NULL){
  
  if (is.null(sample_names)){
    sample_names <- rownames(sampleData(object))
  }  
  
  upstream <- sampleData(object)$upstream
  downstream <- sampleData(object)$downstream
  target_length <- sampleData(object)$body
  bin_size <- sampleData(object)$bin.size
  
  upstreamBins <- upstream/bin_size
  downstreamBins <- downstream/bin_size
  target_bins <- target_length/bin_size
  
  
  enrichMAT <- list()
  for (i in seq_along(assays(object))){
    enrichMAT[[i]] <- assays(object)[[i]]
    
    if(upstreamBins[i] == 0){
      attr(enrichMAT[[i]], "upstream_index") <- vector(mode = "numeric")
    }else{
      attr(enrichMAT[[i]], "upstream_index") <- seq_len(upstreamBins[i])
    }
    
    if(target_bins[i] == 0){
      attr(enrichMAT[[i]], "target_index") <- NULL
      attr(enrichMAT[[i]], "target_is_single_point") <- TRUE
    }else{
      attr(enrichMAT[[i]], "target_index") <- (upstreamBins[i]+1):(upstreamBins[i]+target_bins[i])
      attr(enrichMAT[[i]], "target_is_single_point") <- FALSE
    }
    
    if(downstreamBins[i] == 0){
      attr(enrichMAT[[i]], "downstream_index") <- vector(mode = "numeric")
    }else{
      attr(enrichMAT[[i]], "downstream_index") <- (upstreamBins[i]+target_bins[i]+1):(upstreamBins[i]+target_bins[i]+downstreamBins[i])
    }
    attr(enrichMAT[[i]], "extend") = c(upstream[i], downstream[i]) # it must be a vector of length two
    attr(enrichMAT[[i]], "signal_name") = sample_names[i]
    attr(enrichMAT[[i]], "target_name") = "target"
    class(enrichMAT[[i]]) = c("normalizedMatrix", "matrix")
  }
  names(enrichMAT) <- sample_names
  
  return(enrichMAT)
})


#' generateEnrichedHeatmap
#'
#' export a profileplyr object directly to an object of the EnrichedHeatmap class
#'
#' @rdname generateEnrichedHeatmap
#' @param object A profileplyr object
#' @param include_group_annotation If TRUE (default value) then the Heatmap will be grouped based on the range metadata column specified by 'rowGroupsInUse'
#' @param extra_annotation_columns A character vector of names that match column names of mcols(object). Extra annotation columns will be added to the heatmap based on the values of these indicated range metadata columns.
#' @param sample_names A character vector that will set the names of the heatmap components that are generated from the profileplyr assays() matrices. This argument is optional, by default the names will be the name of the samples in the profileplyr object metadata(proplyrObject)$sampleData$sample_labels.
#' @param return_ht_list Whether the returned object is the heatmap list and not the actual figure. This will be a list of the various components (heatmaps and annotation columns) that can be added to with additional columns in a customized manner.
#' @param ylim A numeric vector of two numbers that species the minimum and maximum of the yaxis of all the heatmaps generated for the matrices. The default is to use the max of the heatmap with the highest signal. If ylim = NULL, different ranges will be inferred for each heatmap. 
#' @param decreasing If object@params$mcolToOrderBy has been changed and is NULL, then the ranges will be ordered by the column indicated in this slot of the metadata. By default, the order will be increasing for the factor or numeric value. For decreasing order, choose decreasing = TRUE.
#' @param all_color_scales_equal If TRUE (default value) then the same color scale will be used for each separate heatmap. If FALSE, color scales will be inferred for each heatmap as indicated by the legends.
#' @param matrices_color Either a single character vector, a numeric vector, a function call to colorRamp2 from the circlize package, or a list. For anything but a list, all the heatmaps generated for the matrices of the profileplyr object will be the same and will be colored as specified here. The character and numeric vector inputs must be either two or three elements in length (denoting color progressions - three elements will give a middle color break), and each element must be a character string or number that points to a color. By default, numeric vectors use the colors in palette(), however this can be expanded with other R color lists(e.g. colors()). If this argument is a list then it's length must equal the number of matrices/samples that exist in the input profileplyr object. The components of the list can be either a numeric vector, character vector, or color function (they do not have to all be the same type of specification). Each element in the list will be the color mapping to the corresponding element in the profileplyr object.
#' @param matrices_pos_line A logical for whether to draw a vertical line(s) at the position of the target (for both a single point or a window). Default is true.
#' @param matrices_pos_line_gp Graphics parameters for the vertical position lines. Should be set with the gpar() function from the grid() package.
#' @param matrices_show_heatmap_legend Logical denoting whether legends for all the heatmaps showing signal over the ranges/matrices should be shown. Default is FALSE.
#' @param matrices_column_title_gp Graphics parameters for the titles on top of each range/matrix (set by 'sample_names' argument or the names of each matrix by default). Should be set with the gpar() function from the grid() package.
#' @param matrices_axis_name_gp Graphics parameters for the text on the x-axis of each matrix heatmap. Should be set with the gpar() function from the grid() package.
#' @param group_anno_color This will specify colors for the grouping column if the 'include_group_annotation' argument is set to TRUE. Since the group column of the range metadata should always be a discrete value, this should be either a numeric vector or character vector with color names. By default, numeric vectors use the colors in palette(), however this can be expanded with other R color lists(e.g. colors()). The length of this vector must equal the number of groups.
#' @param group_anno_width A numeric value that is used to will set the width of the column bar (in mm using the unit() function from the grid package) for the grouping annotation column. 
#' @param group_anno_row_title_gp Graphics parameters for the labels of the groups on the side of the heatmap. Should be set with the gpar() function from the grid() package. 
#' @param group_anno_column_names_gp Graphics parameters for the label of the grouping annotation column. Should be set with the gpar() function from the grid() package.
#' @param extra_anno_color This will specify colors for the annotation columns added by the 'extra_annotation_columns' argument. This must be a list that is of equal length to the 'extra_annotation_columns' argument. Each element of this list will be used to specify the color scheme for the corresponding element of the 'extra_annotation_columns' vector. If an element is NULL, the default colors will be used for the column annotation. For a column with discrete variables this will typically be a vector of numbers or a vector of color names. By default, numeric vectors use the colors in palette(), however this can be expanded with other R color lists(e.g. colors()). For columns with continuous variables, this can also be a a vector of numbers or a vector of color names to signify the color progression, or it can be color mapping function using colorRamp2() from the circlize package.
#' @param extra_anno_top_annotation This is a logical vector that determines whether annotation plots are shown on top of the heatmaps for the extra annotations. This must either be a length of 1, in which case all of the heatamps will abide by this value. Otherwise this must be a vector of equal length to the 'extra_annotation_columns' argument and the elements of this vector will correspond to the equvalent elements in 'extra annotation_columns'
#' @param extra_anno_width This will set the width of the individual extra annotation columns on the right side of the figure. This must be a numeric vector with each element setting the width for the corresponding element in the 'extra_annotation_columns' argument.
#' @param gap The size of the gap between heatmaps and annotations. Only relevant if return_ht_list = FALSE
#' @param genes_to_label A character vector of gene symbols that should match character strings in the 'SYMBOL' column that results from either 'annotateRanges' or 'annotateRanges_great'. Genes that are both in this vector and in the 'SYMBOL' column will be labeled on the heatmap. 
#' @details Takes a profileplyr object and generated heatmap that can be annotated by group or by range metadata columns of the profileplyr object
#' @return By default a customized version of a heatmap from EnrichedHeatmap, if return_ht_list = TRUE then a heatmap list is returned that can be modified and then entered as an input for the \code{\link[EnrichedHeatmap]{EnrichedHeatmap}} function
#' @examples
#' example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
#' object <- import_deepToolsMat(example)
#' 
#' generateEnrichedHeatmap(object, include_group_annotation = FALSE)
#' @export

generateEnrichedHeatmap <- function(object, include_group_annotation = TRUE, extra_annotation_columns = NULL, sample_names = NULL, return_ht_list = FALSE, ylim = "common_max", 
                                    decreasing = FALSE, all_color_scales_equal = TRUE, matrices_color, matrices_pos_line = TRUE, matrices_pos_line_gp = gpar(lty = 2), 
                                    matrices_show_heatmap_legend = TRUE, matrices_column_title_gp =  gpar(fontsize = 10, fontface = "bold"), matrices_axis_name_gp = gpar(fontsize = 8), 
                                    group_anno_color = NULL, group_anno_width = 3, group_anno_row_title_gp = gpar(fontsize = 10), group_anno_column_names_gp = gpar(fontsize = 10),
                                    extra_anno_color = vector(mode = "list", length = length(extra_annotation_columns)), extra_anno_top_annotation = TRUE, 
                                    extra_anno_width = (rep(6, length(extra_annotation_columns))), gap = 2, genes_to_label = NULL){

  # want to order by group, and then within each group order by mean signal
  scoreMat <- do.call(cbind,
                      as.list(assays(object)))
  
  means <- rowMeans(scoreMat)
  means_rev <- max(means) - means
  rowRanges(object)$bin_means_rev <- means_rev
  
  
  rowGroupsInUse <- params(object)$rowGroupsInUse
  mcolToOrderBy <- params(object)$mcolToOrderBy
  if(is.null(mcolToOrderBy)){
    order <- order(mcols(object)[, rowGroupsInUse],
                   mcols(object)$bin_means_rev)
  } else {
    order <- order(mcols(object)[, rowGroupsInUse],
                   mcols(object)[, mcolToOrderBy],
                   decreasing = decreasing)
  }
  
  object <- object[order, ]
  
  enrichMAT <- convertToEnrichedHeatmapMat(object)
  
  # this will find the max for the mean of the columns (bins) and set this as the max for 'ylim' if specified above
  # if ther eis grouping we actualyl want the max ylim necessary when the matrices are divided by group
  
  scoreMat <- do.call(cbind,
                      as.list(enrichMAT))
  
  group_boundaries <- getGroupInfoFromObject(object)$group_boundaries
  # divide the maritcies into separate matrcies for each group
  group_sub <- vector(mode = "list", length = length(group_boundaries)-1)
  for(i in seq_along(group_boundaries[-(length(group_boundaries))])){
    if(i==1){
      group_sub[[i]] <- scoreMat[group_boundaries[i]:group_boundaries[i+1], ]
    }else{
      group_sub[[i]] <- scoreMat[(group_boundaries[i]+1):group_boundaries[i+1], ]
    }
  }
  
  # get the max and min col mean accounting for groups
  col_means <- vector(mode = "list", length = length(group_sub))
  for (i in seq_along(group_sub)){
    col_means[[i]] <- colMeans(group_sub[[i]])
  }
  col_means_unlist <- unlist(col_means)
  
  col_means_max <- max(col_means_unlist)
  # get the value for which this would be 90% to giv esome room in the figure
  max_for_figure = col_means_max/0.9
  
  col_means_min <- min(col_means_unlist)
  min_for_figure <- col_means_min/0.9
  
  if (min_for_figure < 0){
    min_for_figure <- col_means_min/0.9
  }else {
    min_for_figure <- 0
  }
  
  yaxis_side <- vector(length = length(enrichMAT))
  yaxis <- vector(length = length(enrichMAT))
  if(!is.null(ylim)){
    if(ylim == "common_max"){
      ylim <- c(min_for_figure, max_for_figure)
    }
    for (i in seq_along(enrichMAT)){
      if(i == 1){
        yaxis_side[i] <- "left"
        yaxis[i] = TRUE
      } else{
        yaxis_side[i] <- "right"
        yaxis[i] = FALSE
      }
    }
  }else if (is.null(ylim)){
    for (i in seq_along(enrichMAT)){
      yaxis_side[i] <- "right"
      yaxis[i] = TRUE
    }
  }
  
  # following three if statements will show one legend if all color scales are the same and ther user wants a legend (Default)
  # will show all legend if the color scales are unequal and user wants legends

  if(all_color_scales_equal == TRUE & matrices_show_heatmap_legend == TRUE){
    matrices_show_heatmap_legend <- vector(length = length(enrichMAT))
    heatmap_legend_param <- vector(length = length(enrichMAT))
    for (i in seq_along(enrichMAT)){
      if(i == 1){
        matrices_show_heatmap_legend[i] <- TRUE
        heatmap_legend_param[i] <- "signal"
      }
      else{
        matrices_show_heatmap_legend[i] <- FALSE
        heatmap_legend_param[i] <- names(enrichMAT[i])
      }
    }
  } else if(all_color_scales_equal == FALSE & matrices_show_heatmap_legend == TRUE){
    matrices_show_heatmap_legend <- vector(length = length(enrichMAT))
    heatmap_legend_param <- vector(length = length(enrichMAT))
    for (i in seq_along(enrichMAT)){
      matrices_show_heatmap_legend[i] <- TRUE
      heatmap_legend_param[i] <- names(enrichMAT[i])
    }
  } else if(matrices_show_heatmap_legend == FALSE){
    matrices_show_heatmap_legend <- vector(length = length(enrichMAT))
    heatmap_legend_param <- vector(length = length(enrichMAT))
    yaxis_facing <- vector(length = length(enrichMAT))
    for (i in seq_along(enrichMAT)){
      matrices_show_heatmap_legend[i] <- FALSE
      heatmap_legend_param[i] <- names(enrichMAT[i])
    }
  }
  

  if(missing(matrices_color)){
    matrices_color <- vector(mode = "list", length = length(enrichMAT))
  }
  
  # make one big matrix as a means to make common scales across multiple heatmaps.
  # The q_all value will be used below
  all_mats <- unlist(enrichMAT)
  q_all <- quantile(all_mats, 0.99)
  matrices_color_all <- (colorRamp2(c(0, q_all[1]/2, q_all[1]), 
                                     c("blue", "white", "red")))
  
  # this will be run if 'matrices_color' argument is a list of equal length to the number of matrices
  # this will fill in any NULL values in the color list with the default, use a function if its the input, and then put an input vector into the color slot in the default function 
  # if 'all_color_scales_equal' is TRUE, then all heatmaps will be scaled the same based on the quantiles of all values put together (calculated above)
  if (is.list(matrices_color)){
    if(!(length(matrices_color) == length(enrichMAT))){
      stop("The length of the list for the 'matrices_color' argument must equal the number of matrices/samples in the profileplyr object")
    }
    if(all_color_scales_equal == FALSE){
      for(i in seq_along(matrices_color)){
        if(is.null(matrices_color[[i]])){
          q <- quantile(enrichMAT[[i]], 0.99)
          matrices_color[[i]] <- (colorRamp2(c(0, q[1]/2, q[1]), 
                                              c("blue", "white", "red")))
        }else if(is.function(matrices_color[[i]])){
          matrices_color[[i]] <- matrices_color[[i]]
        }else if(is.character(matrices_color[[i]]) | is.numeric(matrices_color[[i]])){
          q <- quantile(enrichMAT[[i]], 0.99)
          if(length(matrices_color[[i]]) == 3){
            matrices_color[[i]] <- colorRamp2(c(0, q[1]/2, q[1]), 
                                               matrices_color[[i]])
          }else if(length(matrices_color[[i]]) == 2){
            matrices_color[[i]] <- colorRamp2(c(0, q[1]), 
                                               matrices_color[[i]])
          }
        }else{
          stop("matrix color specifications must either be a colorRamp2() function call or a character/numeric vector")
        }
      }
    }
    if(all_color_scales_equal == TRUE){
      for(i in seq_along(matrices_color)){
        if(is.null(matrices_color[[i]])){
          matrices_color[[i]] <- (colorRamp2(c(0, q_all[1]/2, q_all[1]), 
                                              c("blue", "white", "red")))
        }else if(is.function(matrices_color[[i]])){
          matrices_color[[i]] <- matrices_color[[i]]
          warning("User provided color function will override the 'all_color_scales = TRUE' option")
        }else if(is.character(matrices_color[[i]]) | is.numeric(matrices_color[[i]])){
          if(length(matrices_color[[i]]) == 3){
            matrices_color[[i]] <- colorRamp2(c(0, q_all[1]/2, q_all[1]), 
                                               matrices_color[[i]])
          }else if(length(matrices_color[[i]]) == 2){
            matrices_color[[i]] <- colorRamp2(c(0, q_all[1]), 
                                               matrices_color[[i]])
          }
        }else{
          stop("matrix color specifications must either be a colorRamp2() function call or a character/numeric vector")
        }
      }
    }
  }
 
  # if the 'matrices_colors' argument is not a list and is a single character or numeric vector, or a function, then all heatmaps will get this scheme, so a list will be populated with the same entries for all
  if (is.character(matrices_color) | is.numeric(matrices_color) | is.function(matrices_color)){
    temp <- vector(mode = "list", length = length(enrichMAT))
    if(all_color_scales_equal == FALSE){
    for(i in seq_along(temp)){
      if(is.function(matrices_color)){
        temp[[i]] <- matrices_color
      }else if(is.character(matrices_color) | is.numeric(matrices_color)){
        q <- quantile(enrichMAT[[i]], 0.99)
        if(length(matrices_color) == 3){
          temp[[i]] <- colorRamp2(c(0, q[1]/2, q[1]),matrices_color)
        }else if(length(matrices_color) == 2){
          temp[[i]] <- colorRamp2(c(0, q[1]), matrices_color)
        }
      }
    }
    }
    if(all_color_scales_equal == TRUE){
      for(i in seq_along(temp)){
        if(is.function(matrices_color)){
          temp[[i]] <- matrices_color
          warning("User provided color function will override the 'all_color_scales = TRUE' option")
        }else if(is.character(matrices_color) | is.numeric(matrices_color)){
          if(length(matrices_color) == 3){
            temp[[i]] <- colorRamp2(c(0, q_all[1]/2, q_all[1]),matrices_color)
          }else if(length(matrices_color) == 2){
            temp[[i]] <- colorRamp2(c(0, q_all[1]), matrices_color)
          }
        }
      }
    }
    matrices_color <- temp
  }
  
  if (is.null(sample_names)){
    sample_names <- names(enrichMAT)
  }  

  if (!(is.null(genes_to_label))){
  # add text row annotation of genes to label
  genes_overlap <- mcols(object)$SYMBOL
  
  index <- vector()
  for (i in seq_along(genes_overlap)){
    if(genes_overlap[i] %in% genes_to_label){
      index <- c(index, i)
    }
  }
  gene_annotation <- rowAnnotation(labels = anno_mark(at = index, 
                                                      labels = genes_overlap[index],
                                                      labels_gp = gpar(fontsize = 6)))
  } else{
    gene_annotation = NULL
  }
  
  heatmap_list <- list()
  
  if (include_group_annotation == TRUE){
    if(is.null(group_anno_color)){
      group_anno_color = seq_along(table((mcols(object)[params(object)$rowGroupsInUse])))
    }
    if(!(length(group_anno_color) == length(seq_along(table((mcols(object)[params(object)$rowGroupsInUse])))))){
      stop("The length of the 'group_anno_color' argument is not the same length as the number of groups in the 'rowGroupsInUse' column of the profileplyr object range metadata.")
    }
    heatmap_list[[1]] <- Heatmap(mcols(object)[,params(object)$rowGroupsInUse], 
                                 col = structure(group_anno_color,
                                                 names = names(table((mcols(object)[params(object)$rowGroupsInUse])))),  
                                 name = params(object)$rowGroupsInUse,
                                 show_row_names = FALSE, 
                                 width = unit(group_anno_width, "mm"),
                                 row_title_gp = group_anno_row_title_gp,
                                 column_names_gp = group_anno_column_names_gp) 

    for (i in (1 + seq_along(enrichMAT))){

      default_color_fun = function(x) {
        q = quantile(x, c(0.01, 0.99))
        qUpper = q[2]
        colorRamp2(c(0, qUpper/2, qUpper), c("blue", "white", "red"))
      }
      
      heatmap_list[[i]] <- EnrichedHeatmap(enrichMAT[[i-1]], 
                                           col = matrices_color[[i-1]],
                                           name = sample_names[i-1],
                                           column_title = sample_names[i-1],
                                           top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = group_anno_color),
                                                                                                    ylim = ylim, 
                                                                                                    axis_param = list(side = yaxis_side[i-1],
                                                                                                                      facing = "outside"),
                                                                                                    yaxis = yaxis[i-1])),
                                           show_heatmap_legend = matrices_show_heatmap_legend[i-1],
                                           heatmap_legend_param = list(title = heatmap_legend_param[i-1]),
                                           column_title_gp = matrices_column_title_gp,
                                           axis_name_gp = matrices_axis_name_gp,
                                           pos_line = matrices_pos_line,
                                           pos_line_gp = matrices_pos_line_gp)
    }
  }
  
  if (include_group_annotation == FALSE){
    for (i in seq_along(enrichMAT)){
      heatmap_list[[i]] <- EnrichedHeatmap(enrichMAT[[i]], 
                                           col = matrices_color[[i]],
                                           name = sample_names[i],
                                           column_title = sample_names[i],
                                           top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = group_anno_color),
                                                                                                    ylim = ylim,
                                                                                                    axis_param = list(side = yaxis_side[i],
                                                                                                                      facing = "outside"),
                                                                                                    yaxis = yaxis[i])),
                                           show_heatmap_legend = matrices_show_heatmap_legend[i],
                                           heatmap_legend_param = list(title = heatmap_legend_param[i]),
                                           column_title_gp = matrices_column_title_gp,
                                           axis_name_gp = matrices_axis_name_gp,
                                           pos_line = matrices_pos_line,
                                           pos_line_gp = matrices_pos_line_gp)
    }
  }
  
  # this will add gene annotation labels to the last heatmap only is there are no additional annotation columns
  if(is.null(extra_annotation_columns) & !(is.null(genes_to_label))){
    if (include_group_annotation == TRUE){
      i = length(enrichMAT) + 1
    }
    if (include_group_annotation == FALSE){
      i = length(enrichMAT)
    }
    
    heatmap_list[[i]] <- EnrichedHeatmap(enrichMAT[[i-1]], 
                                         col = matrices_color[[i-1]],
                                         name = sample_names[i-1],
                                         column_title = sample_names[i-1],
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = group_anno_color),
                                                                                                  ylim = ylim, 
                                                                                                  axis_param = list(side = yaxis_side[i-1],
                                                                                                                    facing = "outside"),
                                                                                                  yaxis = yaxis[i-1])),
                                         show_heatmap_legend = matrices_show_heatmap_legend[i-1],
                                         heatmap_legend_param = list(title = heatmap_legend_param[i-1]),
                                         column_title_gp = matrices_column_title_gp,
                                         axis_name_gp = matrices_axis_name_gp,
                                         pos_line = matrices_pos_line,
                                         pos_line_gp = matrices_pos_line_gp,
                                         right_annotation = gene_annotation)
    
  }
  
  
  if (!is.null(extra_annotation_columns)){
    if (length(extra_anno_top_annotation) == 1){
      extra_anno_top_annotation <- rep(extra_anno_top_annotation, length(extra_annotation_columns))
    } else{
      if (!(length(extra_anno_top_annotation) == length(extra_annotation_columns))){
        stop("Length of 'extra_anno_top_annotation' argument must have a length of one or equal to the length of 'extra_annotation_columns'")
      }
    }
    
    if(!is(extra_anno_color,"list")){
      stop("The 'extra_anno_color' argument must be a list")
    }
    
    if(!(length(extra_anno_color) == length(extra_annotation_columns))){
      stop("The length of the 'extra_annotation_columns' charcter vector must be the same length as the 'extra_anno_color' list")
    }
    
    for (i in seq_along(extra_annotation_columns)){
      
      if(!(extra_annotation_columns[i] %in% colnames(mcols(object)))){
        stop(paste0("'", extra_annotation_columns[i], "' is not a column title in the range metadata."))
      }
      
      column <- mcols(object)[colnames(mcols(object)) %in% extra_annotation_columns[i]][,1] # this is the contents of the column to be visualized, assigning to a variable as this is used often below
      
      if (is.factor(column) | is.character(column)){
        
        if(is.null(extra_anno_color[[i]])){
          extra_anno_color[[i]] = seq_along(table(as.character(column)))
        }
        if(!(length(extra_anno_color[[i]]) == length(seq_along(table(as.character(column)))))){
          stop("The length of the 'extra_anno_color' argument for the '", extra_annotation_columns[i], "' column is not the same length as the number of discrete names in that column of the range metadata.")
        }
        if (extra_anno_top_annotation[i]){
          top_annotation <- HeatmapAnnotation(summary = anno_summary(axis = FALSE))
        }else {
          top_annotation <- NULL
        }
      }else if (is.numeric(column)){
        
        if(is.null(extra_anno_color[[i]])){
          extra_anno_color[[i]] = c("white", "red")
        }
        
        if (extra_anno_top_annotation[i]){
          top_annotation <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = group_anno_color),
                                                                     outline = FALSE, 
                                                                     axis = FALSE))
        } else {
          top_annotation <- NULL
        }

      }else {
        stop(" additional column to be used for annotation must be a numeric, character, or factor variable")
      }
      
      heatmap_list[[length(heatmap_list)+1]] <- Heatmap(column,
                                                        col = extra_anno_color[[i]],
                                                        name = extra_annotation_columns[i],
                                                        show_row_names = FALSE, 
                                                        width = unit(extra_anno_width[i], "mm"),
                                                        column_names_gp = gpar(fontsize = 10),
                                                        top_annotation = top_annotation)
    }
    if(!is.null(genes_to_label)){
      i = length(extra_annotation_columns)
      heatmap_list[[length(heatmap_list)]] <- Heatmap(column,
                                                        col = extra_anno_color[[i]],
                                                        name = extra_annotation_columns[i],
                                                        show_row_names = FALSE, 
                                                        width = unit(extra_anno_width[i], "mm"),
                                                        column_names_gp = gpar(fontsize = 10),
                                                        top_annotation = top_annotation,
                                                        right_annotation = gene_annotation)
    }
  }
  
  
  
  ht_list <- Reduce("+", heatmap_list) 
  
  ifelse(return_ht_list, return(ht_list), return(draw(ht_list, 
                                                      gap =  unit(gap, "mm"),
                                                      split = mcols(object)[,params(object)$rowGroupsInUse])))
  
  #return(draw(ht_list, split = mcols(object)[,object@params$rowGroupsInUse]))
  }


#### From S4Vectors:::selectSome

selectSome <- function(obj, maxToShow = 5, ellipsis = "...", ellipsisPos = c("middle", 
                                                                "end", "start"), quote = FALSE) 
{
  if (is.character(obj) && quote) 
    obj <- sQuote(obj)
  ellipsisPos <- match.arg(ellipsisPos)
  len <- length(obj)
  if (maxToShow < 3) 
    maxToShow <- 3
  if (len > maxToShow) {
    maxToShow <- maxToShow - 1
    if (ellipsisPos == "end") {
      c(head(obj, maxToShow), ellipsis)
    }
    else if (ellipsisPos == "start") {
      c(ellipsis, tail(obj, maxToShow))
    }
    else {
      bot <- ceiling(maxToShow/2)
      top <- len - (maxToShow - bot - 1)
      nms <- obj[c(seq_len(bot), top:len)]
      c(as.character(nms[seq_len(bot)]), ellipsis, as.character(nms[-c(seq_len(bot))]))
    }
  }
  else {
    obj
  }
}






#' Import ChIPprofile object to profileplyr
#'
#' Function to convert soGGi ChIPprofile objects to  profileplyr object .
#'
#' @rdname as_profileplyr
#' @param chipProfile A ChIPprofile object as created by soGGi regionPlot() function.
#' @param names Column to select row IDs/names from ChIPprofile mcols.
#' @return A profileplyr object
#' @examples
#' 
#' library(soGGi)
#' data("ik_Profiles")
#' proplyr <- as_profileplyr(ik_Profiles,names="ID")
#' export_deepToolsMat(proplyr,con=file.path(tempdir(),"ik_Profiles.MAT"))
#' @importClassesFrom soGGi  ChIPprofile
#' @export
as_profileplyr <- function(chipProfile,names = NULL){
  if(!is(chipProfile,"ChIPprofile")) stop("Object must be ChIPprofile object")
  
  if(is.null(names)){
    names <- colnames(mcols(chipProfile))[1]
  }
  
  forDP_Assays <- assays(chipProfile)
  names(forDP_Assays) <-  sample_labels <- basename(metadata(chipProfile)$names)
  
  forDP_ranges <- rowRanges(chipProfile)
  mcols(forDP_ranges) <- cbind(as.data.frame(mcols(forDP_ranges)),names=mcols(forDP_ranges)[,names])
  
  if(!("sgGroup" %in% colnames(mcols(forDP_ranges)))){
    mcols(forDP_ranges)$sgGroup <- "no_groups"
  }
  
  rowGroupsInUse <- list(rowGroupsInUse="sgGroup")
  info <- list(verbose=FALSE,
               scale=1,
               `skip zeros`=FALSE,
               `nan after end`=FALSE,
               `sort using`="mean",
               `unscaled 5 prime`=rep(0,length(forDP_Assays)),
               body=rep(0,length(forDP_Assays)),
               sample_labels=sample_labels,
               downstream=rep(ceiling(ncol(forDP_Assays[[1]])/2),length(forDP_Assays)),
               `unscaled 3 prime`=rep(0,length(forDP_Assays)),
               group_labels=unique(mcols(forDP_ranges)$sgGroup),
               `bin size`=rep(1,length(forDP_Assays)),
               upstream=rep(floor(ncol(forDP_Assays[[1]])/2),length(forDP_Assays)),
               group_boundaries=NA,
               `max threshold`=NULL,
               `ref point`=rep("center",length(forDP_Assays)),
               `min threshold`=NULL,
               `sort regions`="keep",
               `proc number`=11,
               `bin avg type`="mean",
               `missing data as zero`=FALSE)

  
  info_for_sampleData <- info[!(names(info) %in% c("group_labels", "group_boundaries"))] # Doug added this - I think its better to remove right away because group labels doesn't get included if theres multiple groups so it threw an error when you removed it afterwards (removed that code)
  if(length(sample_labels) > 1){
    sampleData <- DataFrame(as.data.frame(info_for_sampleData[lengths(info_for_sampleData) == length(sample_labels)]) %>%
                              cbind(info_for_sampleData %>% .[lengths(.) == 1] %>% lapply(rep,length(sample_labels)) %>% as.data.frame
                              ),
                            row.names = sample_labels)
  }else{
    sampleData <- DataFrame(as.data.frame(info_for_sampleData %>% .[lengths(.) == 1]
    ),
    row.names = sample_labels) 
  }
  sampleData$max.threshold <- sampleData$min.threshold <- replicate(nrow(sampleData),NULL)
  
  
  params <- list(perSampleDPParams= standard_DPparams()$perSampleDPParams,
                 rowGroupsInUse="sgGroup",
                 mcolToOrderBy="sgGroup")
  
  profileplyrDataset <- profileplyr_Dataset(forDP_Assays,forDP_ranges,sampleData,sampleParam=sampleData,params=params)
  return(profileplyrDataset)
}

#' BamBigwig_to_chipProfile
#' 
#' Generate a soGGi ChIPprofile object with multiple BAM/bigWig files or multiple BED files as the input
#'
#'
#' @rdname BamBigwig_to_chipProfile
#' @param signalFiles paths to either BAM files or bigwig files. More than one path can be in this character vector, but all paths in one function call must point to be either all BAM files or all bigWig files, not a combination of the two.
#' @param testRanges Either a character vector with paths to BED files.
#' @param format character vector of "bam", "bigwig", "RleList" or "PWM"
#' @param ... pass to regionPlot() within the soGGi package
#' @return A profileplyr object
#' @examples
#' signalFiles <- c(system.file("extdata",
#'                              "Sorted_Hindbrain_day_12_1_filtered.bam",
#'                               package = "profileplyr"))
#'require(Rsamtools)
#'for (i in seq_along(signalFiles)){
#'  indexBam(signalFiles[i])
#'}
#' testRanges <- system.file("extdata", 
#'                           "newranges_small.bed", 
#'                           package = "profileplyr")
#' BamBigwig_to_chipProfile(signalFiles, 
#'                          testRanges, 
#'                          format = "bam",
#'                          paired=FALSE,
#'                          style="percentOfRegion",
#'                          )
#' @importFrom BiocParallel bplapply
#' @importFrom soGGi regionPlot
#' @importFrom rtracklayer import.bed import.bw
#' @importFrom GenomeInfoDb seqlevelsStyle<- seqlevelsInUse seqlevels
#' @export
#' 
BamBigwig_to_chipProfile <- function(signalFiles, testRanges, format, ...) {
  
  if (missing(format)){
    stop("'format' argument is missing, it must be entered")
  }
  
  if (format == "bigwig"){
    add_chr <- function(genomeCov){
      if (!(any(grepl("chr", names(genomeCov))))) {
        oddChrom <- grepl("GL", names(genomeCov))
        for (i in seq_along(names(genomeCov))){
          if(!oddChrom[i]){
            names(genomeCov)[i] <- paste0("chr", names(genomeCov)[i])
          } else{
            names(genomeCov)[i] <- names(genomeCov)[i]
          }
        }
      }
      return(genomeCov)
    }
    
    genomeCov_list <- lapply(signalFiles, import.bw, as = "RleList")
    genomeCov_list <- lapply(genomeCov_list, add_chr)
    genomeCov_names <- lapply(genomeCov_list, names)
    common_names <- Reduce(intersect, genomeCov_names)
    
    signalFiles_list <- lapply(genomeCov_list, function(x) x[names(x) %in% common_names])
    #names(signalFiles_list) <- signalFiles
    
    group_labels <- vector()
    testRanges_GR <- GRangesList()
    for(i in seq_along(testRanges)){
      temp <- import.bed(testRanges[i])
      seqlevelsStyle(temp) <- "UCSC"
      temp <- temp[seqnames(temp) %in% common_names]
      temp <- temp 
      seqlevels(temp) <- seqlevelsInUse(temp)
      testRanges_GR[[i]] <- temp
      group_labels[i] <- basename(testRanges[i])
    }
    format <- "rlelist"
  }
  
  if(format == "bam"){
    signalFiles_list <- as.list(signalFiles)
    testRanges_GR <- GRangesList()
    group_labels <- vector()
    for(i in seq_along(testRanges)){
      testRanges_GR[[i]] <- import.bed(testRanges[i])
      seqlevelsStyle(testRanges_GR[[i]]) <- "UCSC"
      group_labels[i] <- basename(testRanges[i])
    }
  }
  
  # Construct the group boundaries
  group_boundaries <- c(0)
  for(i in seq_along(testRanges_GR)){
    group_boundaries <- c(group_boundaries, group_boundaries[i] + length(testRanges_GR[[i]]))
  }
  
  testRanges_GR_unlist <- unlist(testRanges_GR)
  names(testRanges_GR_unlist) <- NULL
  
  ChIPprofile_combined <- bplapply(signalFiles_list, regionPlot, testRanges = testRanges_GR_unlist, format = format, ...)
  ChIPprofile_for_proplyr <- do.call(c,ChIPprofile_combined)
  
  # metadata(ChIPprofile_for_proplyr)$group_boundaries <- group_boundaries
  
  if(is.character(signalFiles)){
    metadata(ChIPprofile_for_proplyr)$names <- signalFiles
  }
  
  rowRanges(ChIPprofile_for_proplyr)$sgGroup <- factor(
    rep(group_labels,
        times=diff(group_boundaries) # Doug changed this too, before was 'each=diff(info$group_boundaries)' and this threw a warning that only first elsemnt was used. I think if you have different group sizes this wouldn't work, and 'times' fixes that
    ),
    levels = group_labels
  )
  return(ChIPprofile_for_proplyr)
}



profileplyr_Dataset <-function(matrix,granges,sampleData,sampleParam,params=NULL){
  tempDou <- SummarizedExperiment(matrix,
                                  rowRanges=granges)
  
  # metadata(tempDou)$info <- c(info)
  if(is.null(params)){
    params <- list(perSampleDPParams= standard_DPparams()$perSampleDPParams,
                   rowGroupsInUse=NULL,
                   mcolToOrderBy=NULL
                  )
  }
  
  proplyDataset <- new("profileplyr", tempDou,
              params=params,
              sampleData=sampleData,
              sampleParams=sampleParam)
  return(proplyDataset)
}


standard_DPparams <- function(){
  perSampleDPParams <- c("upstream","downstream","body","bin.size","ref.point","unscaled.5.prime","unscaled.3.prime","sample_labels")
  perComputeDPParams <- c("verbose","bin.avg.type","missing.data.as.zero","scale","skip.zeros","nan.after.end","proc.number","sort.regions","sort.using","min.threshold","max.threshold")
  DPparams <- list(perSampleDPParams=perSampleDPParams,
                   perComputeDPParams=perComputeDPParams)
  return(DPparams)
}



getGroupInfoFromObject <- function(object){
if(!is.null(params(object)$rowGroupsInUse)){
  group_boundaries <- c(which(!duplicated(rowData(object)[params(object)$rowGroupsInUse]))-1,length(object))
  group_labels <- rowData(object)[params(object)$rowGroupsInUse] %>% 
    as.data.frame %>% 
    .[,1] %>%  
    as.vector %>%
    unique 
  if (length(group_labels)  == 1){ 
    group_labels <- list(group_labels) # if you had more than one group, turning it into a list did not create correct output for deepTools metadata, this seems to fix it
  }
}else{
  group_boundaries <- c(0,length(object))
  group_labels <- "no_groups" %>% list
}
  return(list(group_boundaries=group_boundaries,group_labels=group_labels))
}