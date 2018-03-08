#'@title
#'Show method for the IsoDataSet class.
#'@description
#'\code{show} an IsoDataSet object
#'
#'@param object IsoDataSet class object
#'
#'@return Console output of the object
#'
#'@include IsoDataSet-getters.R
#'@exportMethod show
#'@docType methods
#'@name show
#'@rdname IsoDataSet-show
#'@aliases show,IsoDataSet-method
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com}
#'and Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
#'
#'@examples
#'
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'show(myIsoDataSet)
#'
setMethod(f="show", signature=signature(object="IsoDataSet"),
definition=function(object){
    cat("IsoDataSet \n")
    cat("Isoform Counts: \n", sep=" ")
    show(isoCounts(object)[1:min(3, nrow(isoCounts(object))),]);cat("\n") 
    cat("Gene Counts: \n")
    show(geneCounts(object)[1:min(3, nrow(geneCounts(object))),]);cat("\n") 
    cat("Isoform-Gene relationship: \n")
    show(isoGeneRel(object)[1:min(3, nrow(isoGeneRel(object))), ]);cat("\n") 
    cat("Design Matrix: \n")
    show(expData(object)[1:min(3, nrow(expData(object))), ]);cat("\n") 
    cat("Formula: \n")
    cat("\t"); show(designFormula(object));cat("\n") 
    
})
