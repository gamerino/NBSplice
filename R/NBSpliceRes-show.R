#'@title
#'Show method for the NBSpliceRes class.
#'@description
#'\code{show} an NBSpliceRes object
#'
#'@param object NBSpliceRes class object
#'
#'@return Console output of the object
#'
#'@include NBSpliceRes-getters.R
#'@exportMethod show
#'@docType methods
#'@name show
#'@rdname NBSpliceRes-show
#'@aliases show,NBSpliceRes-method
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com}
#'and Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'
#'data(myDSResults, package="NBSplice")
#'
#'show(myDSResults)
#'
setMethod(f="show", signature=signature(object="NBSpliceRes"),
definition=function(object){
    cat("NBSpliceRes data.frame \n")
    cat("Isoform Counts: \n", sep=" ")
    if(nrow(results(object))>0){
    cat("\t"); show(results(object)[1:5,]);cat("\n") 
    }else{
    cat("\t"); show(results(object));cat("\n") 
    }
})
