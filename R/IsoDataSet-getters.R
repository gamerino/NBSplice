#'@title
#'Getters for an IsoDataSet object.
#'@description
#'Obtain information about IsoDataSet slots, according to the given
#'function call.
#'@param object IsoDataSet class object.
#'@param CPM Logical indicating if expression matrix at the isoform level 
#'should be in CPM scale.
#'
#'@return According to the call one of the following objects can be returned
#'\item{matrix}{Isoform expression matrix.}
#'\item{matrix}{Gene expression matrix.}
#'\item{data.frame}{Experiment information.}
#'\item{data.frame}{Isoforms-gene relationships.}
#'\item{formula}{Formula to be used in the GLM fit.}
#'\item{numeric}{Index of low expressed isoforms.}
#'@seealso \code{\link{IsoDataSet}}
#'@name IsoDataSet-getters
NULL
#'@include IsoDataSet.R
#'@exportMethod isoCounts
#'@docType methods
#'@name isoCounts
#'@rdname IsoDataSet-getters
#'@aliases isoCounts-methods
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting isoform expression counts
#'myCounts<-isoCounts(myIsoDataSet)
#'head(myCounts)
setGeneric(name="isoCounts", def=function(object, CPM=TRUE){
    standardGeneric("isoCounts")
})
#'@name isoCounts
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases isoCounts,IsoDataSet-method
setMethod(f="isoCounts", signature="IsoDataSet", 
definition=function(object, CPM=TRUE){
    if(CPM){
        return(cpm(object@counts))
    }else{
        return(object@counts)
    }
})
#'@exportMethod geneCounts
#'@name geneCounts
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases geneCounts-methods
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting gene expression counts
#'myGeneCounts<-geneCounts(myIsoDataSet)
#'head(myGeneCounts)
setGeneric(name="geneCounts", def=function(object){
    standardGeneric("geneCounts")
})
#'@name geneCounts
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases geneCounts,IsoDataSet-method
setMethod(f="geneCounts", signature="IsoDataSet", 
definition=function(object){
    geneIso<-object@isoGeneRel
    genes<-as.character(unique(geneIso[,"gene_id"]))
    df<-do.call(rbind, lapply(seq_along(genes), function(x){
    aux<-geneIso[geneIso[,"gene_id"]==genes[x], "isoform_id"][1]
    return(object@geneCounts[aux,])    
    }))
    rownames(df)<-genes
    colnames(df)<-colnames(object@counts)
    return(df)
})
#'@exportMethod isoGeneRel
#'@name isoGeneRel
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases isoGeneRel-methods
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting isoform-gene relationships
#'myIsoGeneRel<-isoGeneRel(myIsoDataSet)
#'head(myIsoGeneRel)
setGeneric(name="isoGeneRel", def=function(object){
    standardGeneric("isoGeneRel")
})
#'@name isoGeneRel
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases isoGeneRel,IsoDataSet-method
setMethod(f="isoGeneRel", signature="IsoDataSet", 
definition=function(object){
    return(object@isoGeneRel)
})
#'@exportMethod expData
#'@name expData
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases expData-methods
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting the design matrix
#'myDesignMatrix<-expData(myIsoDataSet)
#'myDesignMatrix
setGeneric(name="expData", def=function(object){
    standardGeneric("expData")
})
#'@name expData
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases expData,IsoDataSet-method
setMethod(f="expData", signature="IsoDataSet", 
definition=function(object){
    return(object@colData)
})
#'@exportMethod designFormula
#'@name designFormula
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases designFormula-methods
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting the model formula
#'myFormula<-designFormula(myIsoDataSet)
#'myFormula
setGeneric(name="designFormula", def=function(object){
    standardGeneric("designFormula")
})
#'@name designFormula
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases designFormula,IsoDataSet-method
setMethod(f="designFormula", signature="IsoDataSet", 
definition=function(object){
    return(object@design)
})
#'@exportMethod lowExpIndex
#'@name lowExpIndex
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases lowExpIndex-methods
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Getting index of low expressed isoforms
#'myLowExpIdx<-lowExpIndex(myIsoDataSet)
#'head(myLowExpIdx)
setGeneric(name="lowExpIndex", def=function(object){
    standardGeneric("lowExpIndex")
})
#'@name lowExpIndex
#'@rdname IsoDataSet-getters
#'@inheritParams isoCounts
#'@aliases lowExpIndex,IsoDataSet-method
setMethod(f="lowExpIndex", signature="IsoDataSet", 
definition=function(object){
    return(object@lowExpIndex)
})