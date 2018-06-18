#'@title
#'Print an NBSpliceRes object.
#'@description
#'Generic \code{print} method for NBSpliceRes class and descendants.
#'
#'@param x NBSpliceRes class object.
#'@param ... Included for generic print compatibility.
#'
#'@return Console output of the object.
#'
#'@include NBSpliceRes-show.R
#'@exportMethod print
#'@docType methods
#'@name NBSpliceRes-print
#'@rdname NBSpliceRes-print
#'@aliases print,NBSpliceRes-method
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com},
#'and Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'
#'data(myDSResults, package="NBSplice")
#'
#'print(myDSResults)
#'
setMethod(f="print",signature=signature(x="NBSpliceRes"),
definition = function(x, ...){
    show(x)
}
)
