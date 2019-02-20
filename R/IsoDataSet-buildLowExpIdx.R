#'@title
#'IsoDataSet low expressed isoforms detection.
#'@description
#'\code{buildLowExpIdx} identifies low expressed isoforms and 
#'stores their indexes to facilitate their identification. Low expression is 
#'defined  combining absolute and relative expression thresholds.
#'@param object IsoDataSet class object.
#'@param colName Character indicating the name of the column in the design 
#'matrix to be considered for differential expression analysis.
#'@param ratioThres Numeric indicating the minimum isoform's relative expression
#'value admitted. If one isoform had expression lower than this threshold in at
#'least one sample, thus it will be ignored for further analysis.
#'@param countThres Numeric indicating the isoform's expression threshold. If 
#'one isoform showed a mean expression value lower than this threshold in at
#'least one experimental condition, thus it will be ignored for further
#'analysis.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'@return IsoDataSet object.
#'@include IsoDataSet-print.R
#'@exportMethod buildLowExpIdx
#'@docType methods
#'@name buildLowExpIdx
#'@rdname IsoDataSet-buildLowExpIdx
#'@import BiocParallel
#'@importFrom BiocParallel bpparam
#'@aliases buildLowExpIdx-methods
#'@seealso \code{\link{IsoDataSet}}
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@@bdmg.com.ar}
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Identification of low expressed isoforms
#'myIsoDataSet<-buildLowExpIdx(myIsoDataSet)
setGeneric(name="buildLowExpIdx", def=function(object, colName="condition", 
    ratioThres=0.01, countThres=1, BPPARAM=bpparam()){
    standardGeneric("buildLowExpIdx")
})
#'@name buildLowExpIdx
#'@rdname IsoDataSet-buildLowExpIdx
#'@aliases buildLowExpIdx,IsoDataSet-method
#'@inheritParams buildLowExpIdx
setMethod(f="buildLowExpIdx", signature=signature(object="IsoDataSet"),
    definition=function(object, colName="condition", ratioThres=0.01, 
    countThres=1, BPPARAM=bpparam()){
#
    iso_cm<-isoCounts(object)
    total_counts<-object@geneCounts
    iso_cm<-cbind(iso_cm, total_counts)
    designMatrix<-expData(object)
    if(!any(colnames(designMatrix)==colName)){
        stop(paste("The design matrix should contain a column called", colName,
            sep=" "))
    }
    condLevs<-levels(designMatrix[,colName])
    if(length(condLevs)<2){
        stop(paste("NBSplice needs at least two levels of the ", colName, 
                    "experimental factor", sep=""))
    }
    samplesCols<-lapply(seq_along(condLevs), function(x){
        return(rownames(designMatrix[designMatrix[,colName]==condLevs[x],]))
    })
    # iso_cpm<-round(iso_cm[,do.call(c,samplesCols)])
    # totalC<-round(iso_cm[,paste(do.call(c,samplesCols), "_All", sep="")])
    iso_cpm<-iso_cm[,do.call(c,samplesCols)]
    totalC<-iso_cm[,paste(do.call(c,samplesCols), "_All", sep="")]
    # 
    idxLowRat<-do.call(c, bplapply(seq_len(nrow(iso_cm)), function(x){
        expAndRat<-do.call(rbind,lapply(seq_along(condLevs), function(i){
        dat<-iso_cpm[x,samplesCols[[i]], drop=FALSE]
        ratios<-dat/totalC[x,paste(samplesCols[[i]], "_All", sep="")]
        return(cbind(meanCond=rowMeans(dat), ratCond=min(ratios)))
    }))
            dat<-iso_cpm[x,]
            allCounts<-totalC[x, , drop=FALSE]
            if(all(dat == allCounts)){ #only one isoform
                # if(rowMeans(allCounts) ==0 | any(expAndRat[,"meanCond"] < 
                #     countThres) ){
                        j<-x
                        # }else{
                        # j<-NULL}
            }else{
                if( any(is.na(expAndRat[,"ratCond"])) | 
                    any(expAndRat[,"meanCond"] < countThres)  ){
                        j<-x
                }else{
                    if(min(expAndRat[,"ratCond"]) < ratioThres){
                        j<-x
                    }else{
                        j<-NULL
                    }
                }
            }

        return(j)
    }, BPPARAM=BPPARAM))
    if(is.null(idxLowRat)) idxLowRat<-numeric()
    object@lowExpIndex<-idxLowRat
    return(object)
})
