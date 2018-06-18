#'@title
#'NBSpliceRes S4 class implementation in R
#'@description 
#'This S4 class represents a data set containing NBTest results results. 
#'
#'@slot results Data.frame with NBTest results of expressed isoforms.
#'@slot lowExpIndex Numeric indicating the positions of low expressed isoforms.
#'@slot contrast Character indicating the contrast used for NBTest.
#'@slot dispersion Numeric model dispersions.
#'@section Features: 
#'\enumerate{
#'  \item Explore differential splicing occurrence. 
#'  \item Explore isoform relative expression an its changes.
#'  \item Combine the results of both gene and isoform levels analysis.
#'}
#'
#'@section Functions:
#'NBSpliceRes S4 class includes the following functions:
#'\describe{
#' \item{initialize}{Constructor of NBSpliceRes objects.} 
#' \item{results}{Gets the results slot.}
#' \item{contrast}{Gets the contrast slot.}
#' \item{lowExpIndex}{Gets the lowExpIdx slot.}
#' \item{disp}{Gets the dispersion slot.}
#' \item{print}{Shows a NBSpliceRes object.}
#' \item{show}{Shows a NBSpliceRes object.}
#' \item{GetDSGenes}{Returns the list of differentially spliced genes.}
#' \item{GetDSResults}{Returns the differential splicing results.}
#' \item{GetGeneResults}{Returns the NBSplice results for an specific gene.}
#' }
#'
#'@include IsoDataSet-core.R
#'@name NBSpliceRes-class
#'@rdname NBSpliceRes-class
#'@exportClass NBSpliceRes
#'@import methods
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#' Fernandez \email{efernandez@@bdmg.com.ar}
#'@examples
#'
#'data(myDSResults, package="NBSplice")
#'myResults<-results(myDSResults, filter=FALSE)
#'myLowExpIdx<-lowExpIndex(myDSResults)
#'myContrast<-contrast(myDSResults)
#'
#'myNewDSResults<-NBSpliceRes(myResults, myLowExpIdx, myContrast)
#'
#'##Getting differentially spliced genes
#'myDSGenes<-GetDSGenes(myDSResults)
#'
#'##Getting the results for differentially spliced genes
#'myDSResultsDF<-GetDSResults(myDSResults)
#'
#'##Getting the results for a particular differentially spliced gene
#'myResults<-results(myDSResults)
#'
#'## Select the first gene
#'gene<-myResults[,"gene"][1]
#'
#'myGeneResults<-GetGeneResults(myDSResults, gene)
#'
setClass(Class="NBSpliceRes", slots=list(results ="data.frame", 
    lowExpIndex="numeric", contrast="character", dispersion="numeric"), 
    validity=function(object){
    ## Check contrast
    contr<-object@contrast

    if(!is.character(contr)){
        stop("The contrast object should be a character")
    }

    ## Check results matrix
    resDF<-results(object, filter=FALSE)
    if(!is.data.frame(resDF)){
        stop(paste("The results object should be a data.frame"))
    }
    if(ncol(resDF)>0 & !all(c("iso", "gene", "odd", "stat", "pval", "genePval",
    "FDR", "geneFDR", paste("ratio", contr, sep="_")) %in% colnames(resDF))){
            stop(paste("The results data.frame should contain at least the ",
                "next column names:","iso,", "gene,", "odd,", "stat,", "pval,",
                "genePval,", "FDR,", "geneFDR,", paste("ratio", contr, 
                sep="_", collapse=" and "), sep=" "))
    }
    ## Check idxLowExp vector
    lowExpIdx<-lowExpIndex(object)
    if(!is.numeric(lowExpIdx)){
        stop("The lowExpIndex object should be a numeric")
    }
    ## Check dispersion vector
    dispersion<-disp(object)
    if(!is.numeric(dispersion)){
        stop("The dispersion object should be a numeric")
    }

    }, prototype=list(
    results =as.data.frame(matrix(ncol=0, nrow=0)),
    lowExpIndex=numeric(),
    contrast=character(),
    dispersion=numeric())
)