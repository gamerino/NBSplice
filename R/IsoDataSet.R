#'@title
#'IsoDataSet S4 class implementation in R
#'@description 
#'This S4 class represents a data set containing isoforms expresion in R. 
#'
#'@slot counts Matrix containing expression values at the isoform level 
#'@slot geneCounts Matrix with expression values at the gene level
#'@slot colData Data.frame with experiment information. Its rows should be 
#' the columns of the counts data.frame
#'@slot isoGeneRel Data.frame specifying the isoform-gene relationship
#'@slot design Formula to be used in the GLM fit and tests
#'@slot lowExpIndex Numeric indicating the positions of low expressed isoforms
#'@section Features: 
#'\enumerate{
#'  \item Discover differential modifications in the splicing patterns. 
#'  \item Detect and quantify isoform relative expression changes.
#'  \item Combine the results of both gene and isoform levels analysis.
#'}
#'@section Functions:
#'IsoDataSet S4 class includes the following functions:
#'\describe{
#' \item{initialize}{Constructor of IsoDataSet objects.} 
#' \item{isoCounts}{Get the counts slot.}
#' \item{geneCounts}{Get the geneCounts slot.}
#' \item{isoGeneRel}{Get the isoGeneRel slot.}
#' \item{colData}{Get the colData slot.}
#' \item{designFormula}{Get the design slot.}
#' \item{lowExpIndex}{Get the lowExpIdx slot.}
#' \item{setDesign}{Set the design slot.}
#' \item{buildLowExpIdx}{Build the index to identify low expressed isoforms.}
#' \item{NBTest}{Perform differential expression analysis.}
#'}
#'@include NBSplice-package.R
#'@name IsoDataSet-class
#'@rdname IsoDataSet-class
#'@exportClass IsoDataSet
#'@import methods
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#' Fernandez \email{efernandez@@bdmg.com.ar}
#'@examples
#'## Data loading
#'data(isoCounts, package="NBSplice")
#'data(designMatrix, package="NBSplice")
#'data(geneIso, package="NBSplice")
#'
#'## Arguments definition
#'colName<-"condition"
#'
#'## Constructor calling
#'myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)
#'
#'## Identificating Low expressed Isoforms
#'myIsoDataSet<-buildLowExpIdx(myIsoDataSet)
#'
#'##Arguments definition
#'colName<-"condition"
#'test<-"F"
#'
#'## Differential splicing test
#'myDSResults<-NBTest(myIsoDataSet, colName, test)
setClass(Class="IsoDataSet", slots=list(counts ="matrix", geneCounts="matrix",
    colData="data.frame", isoGeneRel="data.frame", design="formula", 
    lowExpIndex="numeric"), validity=function(object){

    ## Check counts matrix
    expMat<-isoCounts(object)
    if(!is.matrix(expMat)){
        stop(paste("The argument 'counts' should be a matrix"))
    }
    ## Check geneCounts matrix
    geneExpMat<-object@geneCounts
    if(!is.matrix(geneExpMat)){
        stop("The argument 'geneCounts' should be a matrix")
    }
    ## Check colData matrix
    colDataM<-expData(object)
    if(!is.data.frame(colDataM)){
        stop("The argument 'colData' should be a data.frame")
    }
    if(any(colnames(expMat) != rownames(colDataM))){
        stop("The column names of counts and the row names of geneCounts must be 
            the same")
    }
    
}, prototype=list(
    counts =matrix(ncol=0, nrow=0),
    geneCounts=matrix(ncol=0, nrow=0), 
    colData=data.frame(),
    isoGeneRel=data.frame(),
    design=formula(),
    lowExpIndex=numeric())
)
