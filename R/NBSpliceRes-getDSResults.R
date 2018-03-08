#'@title
#'Get differential expression results of significant genes. 
#'@description
#'\code{GetDSResults} returns the results obtained for those genes identified 
#'as differentially spliced.
#'
#'@param myNBRes NBSpliceRes class object.
#'@param adjusted Logical indicating if adjusted p values should be used.
#'@param p.value Numeric value between 0 and 1 giving the required family-wise
#'error rate or false discovery rate.
#'
#'@return data.frame with the results obtainted by means of the NBTest method
#'
#'@include NBSpliceRes-getDSGenes.R
#'@exportMethod GetDSResults
#'@docType methods
#'@name GetDSResults
#'@rdname NBSpliceRes-GetDSResults
#'@aliases GetDSResults-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'myDSResultsDF<-GetDSResults(myDSResults)
#'
setGeneric(name="GetDSResults", def=function(myNBRes, adjusted=TRUE, 
    p.value=0.05){
        standardGeneric("GetDSResults")
})
#'@name GetDSResults
#'@rdname NBSpliceRes-GetDSResults
#'@aliases GetDSResults,NBSpliceRes-method
#'@inheritParams GetDSResults
setMethod(f="GetDSResults", signature="NBSpliceRes", definition=function(
    myNBRes, adjusted=TRUE, 
    p.value=0.05){
        sigRes<-results(myNBRes, filter=TRUE)
        if(adjusted){
            DSRes<-sigRes[sigRes[, "geneFDR"] < p.value & !is.na(sigRes[, 
                "geneFDR"] ), ]    
            DSRes<-sigRes[sigRes[, "genePval"] < p.value & !is.na(sigRes[, 
            "genePval"]), ]    
        }
        return(DSRes)
})