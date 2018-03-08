#'@title
#'Getters for an NBSpliceRes object.
#'@description
#'Obtain information about NBSpliceRes slots, according to the given
#'function call.
#'@param object NBSpliceRes class object.
#'@param filter Logical indicating if low expressed isoforms should be filtered
#'out.
#'@return according to the call one of the following objects can be returned
#' \item{data.frame}{NBTest results.}
#' \item{numeric}{Index of low expressed isoforms.}
#' \item{character}{Experiment contrast.}
#' \item{numeric}{Estimated gene dispersion.}
#'@name NBSpliceRes-getters
NULL
#'@include NBSpliceRes-constructor.R
#'@exportMethod results
#'@docType methods
#'@name results
#'@rdname NBSpliceRes-getters
#'@aliases results-methods
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'## Getting results slot filtering the low expressed isoforms
#'myResults<-results(myDSResults)
#'
#'## Getting results slot keeping the low expressed isoforms
#'myResults<-results(myDSResults, filter=FALSE)
setGeneric(name="results", def=function(object, filter=TRUE){
    standardGeneric("results")
})
#'@name results
#'@rdname NBSpliceRes-getters
#'@inheritParams results
#'@aliases results,IsoDataSet-method
setMethod(f="results", signature="NBSpliceRes", 
definition=function(object, filter=TRUE){
    if(filter & nrow(object@results)>0){
        ret<-object@results[-object@lowExpIndex,]
        ret[,"pval"]<-as.numeric(as.character(ret[,"pval"]))
        return(ret[!is.na(ret[,"pval"]),])
    }else{
        return(object@results)
    }
})
#'@exportMethod contrast
#'@name contrast
#'@rdname NBSpliceRes-getters
#'@inheritParams results
#'@aliases contrast-methods
#'@examples
#'## Getting the contrast slot
#'myContrast<-contrast(myDSResults)
setGeneric(name="contrast", def=function(object){
    standardGeneric("contrast")
})
#'@name contrast
#'@rdname NBSpliceRes-getters
#'@inheritParams contrast
#'@aliases contrast,NBSpliceRes-method
setMethod(f="contrast", signature="NBSpliceRes", 
definition=function(object){
    return(object@contrast)
})
#'@name lowExpIndex
#'@rdname NBSpliceRes-getters
#'@inheritParams results
#'@aliases lowExpIndex,NBSpliceRes-method
#'@examples
#'## Getting the lowExpIndex slot
#'myLowExpIndex<-lowExpIndex(myDSResults)
#'
setMethod(f="lowExpIndex", signature="NBSpliceRes", 
definition=function(object){
    return(object@lowExpIndex)
})
#'@exportMethod disp
#'@name disp
#'@rdname NBSpliceRes-getters
#'@inheritParams results
#'@aliases disp-methods
#'@examples
#'## Getting the dispersion slot
#'myDispersion<-disp(myDSResults)
setGeneric(name="disp", def=function(object){
    standardGeneric("disp")
})
#'@name disp
#'@rdname NBSpliceRes-getters
#'@inheritParams disp
#'@aliases disp,NBSpliceRes-method
setMethod(f="disp", signature="NBSpliceRes", 
definition=function(object){
    return(object@dispersion)
})