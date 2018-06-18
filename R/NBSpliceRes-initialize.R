#'@title
#'NBSpliceRes object constructor.
#'@description
#'\code{initialize} creates an NBSpliceRes object
#' 
#'@param .Object NBSpliceRes class object.
#'@param results Data.frame with NBTest results of expressed isoforms.
#'@param lowExpIndex Numeric indicating the positions of low expressed isoforms.
#'@param contrast Character indicating the contrast used for NBTest.
#'@param dispersion Numeric with the estimated gene dispersions.
#'
#'@return NBSpliceRes object.
#'
#'@include NBSpliceRes.R
#'@exportMethod initialize
#'@docType methods
#'@name NBSpliceRes-initialize
#'@rdname NBSpliceRes-initialize
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bplapply
#'@aliases initialize,NBSpliceRes-method
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'
#'data(myDSResults, package="NBSplice")
#'myResults<-results(myDSResults, filter=FALSE)
#'myLowExpIdx<-lowExpIndex(myDSResults)
#'myContrast<-contrast(myDSResults)
#'
#'myNewDSResults<-NBSpliceRes(myResults, myLowExpIdx, myContrast)
#'
setMethod(f="initialize", signature=signature(.Object="NBSpliceRes"),
definition=function(.Object, results, lowExpIndex, contrast, dispersion){
    ##Set the different slots
    if(nargs() >=3){
    # results
        .Object@results<-results
        .Object@lowExpIndex<-lowExpIndex
        if(ncol(results)>0 & !(all(paste("ratio", contrast, sep="_") %in% 
            colnames(results)))){
                stop("The result name should be columns called ratio_X, where X
                    means each of the contrasted conditions specified with the
                    contrast parameter.")
        }
        .Object@contrast<-contrast
        if(is.null(names(dispersion))){
            names(dispersion)<-unique(results[,"gene"])
        }
        .Object@dispersion<-dispersion
    }else{

        .Object@results<-as.data.frame(matrix(ncol=0, nrow=0))
        .Object@contrast<-character()
        .Object@lowExpIndex<-numeric()
        .Object@dispersion<-numeric()
    }    
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
    
})
