#'@title
#'IsoDataSet constructor
#'@description
#'\code{IsoDataSet} creates an object to store expression counts at
#'isoform and gene level, the relationship between those, the experiment data
#'and the formula to be used for models fitting required to evaluate 
#'differential splicing.
#'
#'@param isoCounts Matrix having the expression counts at the isoform level. 
#'Isoforms must be in rows and samples in columns. Rownames and colnames must 
#'be defined with isoform and samples names, respectively.
#'@param experimentData Data.frame specifying metadata related to the
#'experiment. Its rows must be the samples and experimental factors should be
#'arranged on its columns.
#'@param colName Character indicating the name of the column in the design 
#'matrix to be considered for differential splicing analysis.
#'@param geneIso Data.frame containing the relationship between isoforms and 
#'genes. It must contain two columns, named as 'gene_id' and 'isoform_id'.
#'Its isoforms should be the same specified in the isoCounts matrix.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return IsoDataSet object.
#'
#'@include IsoDataSet-initialize.R
#'@export IsoDataSet
#'@docType methods
#'@name IsoDataSet
#'@import BiocParallel
#'@importFrom BiocParallel bpparam
#'@rdname IsoDataSet-constructor
#'@aliases IsoDataSet-methods
#'@seealso \code{\link{IsoDataSet-class}}
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com}, and Elmer A. 
#'Fernandez \email{efernandez@@bdmg.com.ar}
#'
#'@examples
#'
#'## Data loading
#'data(isoCounts, package="NBSplice")
#'data(designMatrix, package="NBSplice")
#'data(geneIso, package="NBSplice")
#'
#'## Arguments definition
#'colName<-"condition"
#'
#'## Constructor calling
#'myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)
IsoDataSet<-function(isoCounts, experimentData, colName, geneIso, 
BPPARAM=bpparam()){
    .Object<-new("IsoDataSet")
    if(nargs() >=3 ){
        .Object<-initialize(.Object, isoCounts, experimentData, colName, 
            geneIso, BPPARAM)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
}