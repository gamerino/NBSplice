#'@title
#'Get differentially spliced genes.
#'@description
#'\code{GetGeneResults} returns the results obtained for the specified gene.
#'
#'@param myNBRes NBSpliceRes class object.
#'@param gene Character indicating the gene name. 
#'@param filterLowExpIso Logical indicating if lower-expression isoforms should
#'be filtered out.
#'
#'@return Data.frame object with gene results.
#'
#'@include NBSpliceRes-getDSResults.R
#'@exportMethod GetGeneResults
#'@docType methods
#'@name GetGeneResults
#'@rdname NBSpliceRes-GetGeneResults
#'@aliases GetGeneResults-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'myResults<-results(myDSResults)
#'
#'## Select the first gene
#'gene<-myResults[,"gene"][1]
#'
#'myGeneResults<-GetGeneResults(myDSResults, gene)
#'
setGeneric(name="GetGeneResults", def=function(myNBRes, gene, 
    filterLowExpIso=TRUE){
        standardGeneric("GetGeneResults")
})
#'@name GetGeneResults
#'@rdname NBSpliceRes-GetGeneResults
#'@aliases GetGeneResults,NBSpliceRes-method
#'@inheritParams GetGeneResults
setMethod(f="GetGeneResults", signature="NBSpliceRes", 
    definition=function(myNBRes, gene, filterLowExpIso){
        if(!is.logical(filterLowExpIso)){
            stop("The argument 'filterLowExpIso' should be TRUE or FALSE")
        }
        if(!is.character(gene)){
            stop("The argument 'gene' should be a character")
        }
        sigRes<-results(myNBRes, filter=filterLowExpIso)
        
        if(!any( gene == sigRes[, "gene"]) & !filterLowExpIso){
            stop(paste("The ", gene, " gene is not present in NBTest 
                results.", sep=""))
        }
        
        if(!any( gene == sigRes[, "gene"]) & filterLowExpIso){
            stop(cat(paste("The ", gene, " gene is not present in NBTest 
                results.", sep=""), "\n", "Please, try setting the 
                'filterLowExpIso' to 'FALSE'"))
        }
        if(filterLowExpIso){
            sigRes<-sigRes[!is.na(sigRes[, "pval"]),]
        }
        return(sigRes[sigRes[, "gene"]== gene,])
})