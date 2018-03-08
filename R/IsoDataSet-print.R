#'@title
#'Print an IsoDataSet object.
#'@description
#'Generic \code{print} method for IsoDataSet class and descendants.
#'
#'@param x IsoDataSet class object.
#'@param ... Included for generic print compatibility.
#'
#'@return Console output of the object.
#'
#'@include IsoDataSet-show.R
#'@exportMethod print
#'@docType methods
#'@name print
#'@rdname IsoDataSet-print
#'@aliases print,IsoDataSet-method
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com},
#'and Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
#'
#'@examples
#'
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'print(myIsoDataSet)
#'
setMethod(f="print",signature=signature(x="IsoDataSet"),
definition = function(x, ...){
    show(x)
}
)
