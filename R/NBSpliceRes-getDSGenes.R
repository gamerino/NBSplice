#'@title
#'Get differentially spliced genes.
#'@description
#'\code{GetDSGenes} returns the list of genes identified as differentially 
#'spliced.
#'@param myNBRes NBSpliceRes class object.
#'@param adjusted Logical indicating if adjusted p values should be used.
#'@param p.value Numeric value between 0 and 1 giving the required family-wise
#'error rate or false discovery rate.
#'@return A character with the names of differentially spliced genes.
#'@include IsoDataSet-NBTest.R
#'@exportMethod GetDSGenes
#'@docType methods
#'@name GetDSGenes
#'@rdname NBSpliceRes-GetDSGenes
#'@aliases GetDSGenes-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'myDSGenes<-GetDSGenes(myDSResults)
setGeneric(name="GetDSGenes", def=function(myNBRes, adjusted=TRUE,p.value=0.05){
        standardGeneric("GetDSGenes")
})
#'@name GetDSGenes
#'@rdname NBSpliceRes-GetDSGenes
#'@aliases GetDSGenes,NBSpliceRes-method
#'@inheritParams GetDSGenes
setMethod(f="GetDSGenes", signature="NBSpliceRes", 
    definition=function(myNBRes, adjusted=TRUE, p.value=0.05){
        if(!is.logical(adjusted)){
            stop("The parameter 'adjusted' should be TRUE or FALSE")
        }
        if(!is.numeric(p.value) | p.value < 0 | p.value >1 ){
            stop("The parameter 'p.value' should be a number between 0 and 1")
        }
        
        sigRes<-results(myNBRes, filter=TRUE)
        if(adjusted){
            DSGenes<-sigRes[sigRes[, "geneFDR"] < p.value & !is.na(sigRes[, 
                "geneFDR"] ), "gene"]    
        }else{
            DSGenes<-sigRes[sigRes[, "genePval"] < p.value & !is.na(sigRes[, 
            "genePval"]), "gene"]    
        }
        return(unique(as.character(DSGenes)))
})