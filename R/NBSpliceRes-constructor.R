#'@title
#'NBSpliceRes constructor
#'@description
#'\code{NBSpliceRes} creates an object to store the results of differential
#'expression analysis performed by the NBTest function.
#'
#'@param results Data.frame with NBTest results of expressed isoforms.
#'@param lowExpIndex Numeric indicating the positions of low expressed isoforms.
#'@param contrast Character indicating the contrast used for NBTest.
#'@param dispersion Numeric with the estimated gene dispersions.
#'
#'@return NBSpliceRes object.
#'
#'@include NBSpliceRes-initialize.R
#'@export NBSpliceRes
#'@docType methods
#'@name NBSpliceRes
#'@rdname NBSpliceRes-constructor
#'@aliases NBSpliceRes-methods
#'@seealso \code{\link{NBSpliceRes-class}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com}, and Elmer A. 
#'Fernandez \email{efernandez@@bdmg.com.ar}
#'@examples
#'
#'data(myDSResults, package="NBSplice")
#'myResults<-results(myDSResults, filter=FALSE)
#'myLowExpIdx<-lowExpIndex(myDSResults)
#'myContrast<-contrast(myDSResults)
#'myDispersion<-disp(myDSResults)
#'myNewDSResults<-NBSpliceRes(myResults, myLowExpIdx, myContrast, myDispersion)
#'
NBSpliceRes<-function(results, lowExpIndex, contrast, dispersion){
    .Object<-new("NBSpliceRes")
    if(nargs() >=4 ){
        .Object<-initialize(.Object, results, lowExpIndex, contrast, dispersion)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
}