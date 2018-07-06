#'@title
#'IsoDataSet object constructor.
#'@description
#'\code{initialize} creates an IsoDataSet object with the specified isoform
#'expression counts.
#'
#'@param .Object IsoDataSet class object.
#'@param isoCounts Matrix having the expression counts at the isoform level. 
#'Isoforms must be in rows and samples in columns. Rownames and colnames must 
#'be defined with isoform and samples names, respectively.
#'@param experimentData Data.frame specifying metadata related to the
#'experiment. Its rows must be the samples and experimental factors should be
#'arranged on its columns.
#'@param colName Character indicating the name of the column in the design 
#'matrix to be considered for mean expression calculations per experimental
#'condition and differential expression test.
#'@param geneIso Data.frame containing the relationship between isoforms and 
#'genes. It must contain two columns, named as 'gene_id' and 'isoform_id'.
#'Its isoforms should be the same specified in the isoCounts matrix.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return IsoDataSet object.
#'
#'@include totalGeneCounts.R
#'@exportMethod initialize
#'@docType methods
#'@name initialize
#'@rdname IsoDataSet-initialize
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bplapply
#'@aliases initialize,IsoDataSet-method
#'@seealso \code{\link{IsoDataSet}}
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
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
#'
setMethod(f="initialize", signature="IsoDataSet",
definition=function(.Object, isoCounts, experimentData, colName, geneIso, 
BPPARAM=bpparam()){
    ##Set the different slots
    if(nargs() >=3){
    ## Check isoGeneRel
    if(!any(colnames(geneIso) =="isoform_id")){
        stop("The geneIso data.frame must be two columns called 'gene_id' and
            'isoform_id'. ")
    }
    if(!all(rownames(isoCounts) %in% geneIso[,"isoform_id"])){
        stop("The geneIso object and the isoCounts object should contain the
            same isoforms")
    }
    # counts
        if(is.data.frame(isoCounts)){
            isoCounts<-as.matrix(isoCounts)
            # for(i in seq_len(ncol(isoCounts))){
            #     isoCounts[,i]<-as.numeric(as.character(isoCounts[,i]))
            # }
            isoCounts[,seq_len(ncol(isoCounts))]<-vapply(
                seq_len(ncol(isoCounts)), function(i){
                    return(as.numeric(as.character(isoCounts[,i])))
                }, as.numeric(as.character(isoCounts[,1])))
        }
        if(length(rownames(isoCounts)) ==0 | length(colnames(isoCounts)) ==0) {
            stop("The isoCounts object should have row and column names")
        }
        if(is.factor(geneIso[,"gene_id"])){
            geneIso[,"gene_id"]<-as.character(geneIso[,"gene_id"])
        }
        if(is.factor(geneIso[,"isoform_id"])){
            geneIso[,"isoform_id"]<-as.character(geneIso[,"isoform_id"])
        }
        #sort the IsoCounts and GeneIso by gene name
        geneIso<-geneIso[order(geneIso[,"gene_id"]),]
        isoCounts<-isoCounts[geneIso[,"isoform_id"],]
        .Object@counts<-isoCounts
        # isoGeneRel slot
        .Object@isoGeneRel<-geneIso
        # geneCounts slot
        geneCounts<-totalGeneCounts(isoCounts, geneIso, BPPARAM)
        .Object@geneCounts<-geneCounts
        # colData slot
        colData<-experimentData[colnames(isoCounts), ]
        .Object@colData<-colData
        # design slot
        .Object@design<-as.formula(paste("counts~", colName, "+iso+", colName, 
            ":iso", sep=""))
        .Object@lowExpIndex<-integer()
    }else{

        .Object@counts<-matrix(ncol=0, nrow=0)
        .Object@geneCounts<-matrix(ncol=0, nrow=0)
        .Object@colData<-data.frame()
        .Object@isoGeneRel<-data.frame()
        .Object@design<-formula()
        .Object@lowExpIndex<-integer()
    }    
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
})
