#'@title
#'Differential splicing analysis.
#'@description
#'\code{NBTest} is the main function of NBSplice. It fits the negative binomial
#'models for all genes and performs the corresponding hypothesis tests to
#'evaluate the occurrence of differential splicing. 
#'@param object IsoDataSet class object.
#'@param colName Character indicating the name of the column in the design 
#'matrix to be considered for mean expression calculations per experimental
#'condition.
#'@param test Character indicating the name of the distribution to assume for
#'linear hypothesis statistic. Could be "F" or "chisq".
#'@param contrast Character vector with the names of the two levels of the
#'experimental factor to be contrasted.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'@return A NBSpliceRes object with the obtained results.
#'@include NBSpliceRes-print.R
#'@exportMethod NBTest
#'@docType methods
#'@name NBTest
#'@rdname IsoDataSet-NBTest
#'@import BiocParallel
#'@importFrom BiocParallel bpparam
#'@import stats
#'@importFrom stats p.adjust
#'@aliases NBTest-methods
#'@seealso \code{\link{IsoDataSet}}
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'## Data loading
#'data(myIsoDataSet, package="NBSplice")
#'
#'## Identificating Low expressed Isoforms
#'myIsoDataSet<-buildLowExpIdx(myIsoDataSet)
#'
#'##Arguments definition
#'colName<-"condition"
#'test<-"F"
#'
#'## Differential splicing test
#'myDSResults<-NBTest(myIsoDataSet, colName, test)
#'
setGeneric(name="NBTest", def=function(object, colName="condition", 
    test=c("F", "Chisq"), contrast=c(levels(expData(object)[,colName])[
        seq_len(2)]), BPPARAM=bpparam()){
    standardGeneric("NBTest")
})
#'@name NBTest
#'@rdname IsoDataSet-NBTest
#'@aliases NBTest,IsoDataSet-method
#'@inheritParams NBTest
setMethod(f="NBTest", signature=signature(object="IsoDataSet"),
    definition=function(object, colName="condition", test=c("F", "Chisq"), 
    contrast=c(levels(expData(object)[,colName])[seq_len(2)]), 
    BPPARAM=bpparam()){
    if(!test %in% c("F", "Chisq")){
        stop("The 'test' parameter should be 'F' or 'Chisq'. ")
    }
    if(!colName %in% colnames(expData(object)) | ! grepl(colName, 
        as.character(designFormula(object))[3])){
        stop("The 'colName' parameter should be one of the design matrix 
            columns and the same used to built the object.")
    }
        
    idxLowRat<-lowExpIndex(object)
    # if(length(idxLowRat)==0){
    #     object<-buildLowExpIdx(object)
    #     idxLowRat<-lowExpIndex(object)
    # }
    if(length(idxLowRat)==0){
        warning("The low-expressed isoforms have not been previously
            identified. Please consider their identification before executing
                the NBTest function")
    }
    designMatrix<-expData(object)
    if(length(idxLowRat)>0){
        geneIso<-isoGeneRel(object)[-idxLowRat,, drop=FALSE]
        geneCounts<-object@geneCounts[-idxLowRat,, drop=FALSE]
        isoCounts<-isoCounts(object)[-idxLowRat,, drop=FALSE]
    }else{
        geneIso<-isoGeneRel(object)
        geneCounts<-object@geneCounts
        isoCounts<-isoCounts(object)
    }
    genes<-as.character(unique(geneIso[,"gene_id"]))

    resultsOK<-do.call(rbind, bplapply(seq_along(genes), function(j){
        gene<-genes[j]
        formula<-designFormula(object)
        myData<-buildData(isoCounts, geneCounts, geneIso, gene=gene, 
            designMatrix, colName=colName)
        isoRes<-fitModel(myData, gene=gene,formula=formula, colName=colName,
            test=test, contrast=contrast)
        return(isoRes)
    }, BPPARAM=BPPARAM))
    resultsOK<-as.data.frame(resultsOK)
    # resultsOK<-resultsOK[!is.na(resultsOK$pval),]
    # for(i in seq(from=3, to=ncol(resultsOK))){
    #     resultsOK[,i]<-as.numeric(as.character(resultsOK[,i]))
    # }
    resultsOK[,3:ncol(resultsOK)]<-vapply(3:ncol(resultsOK), function(i){
            return(as.numeric(as.character(resultsOK[,i])))
        }, as.numeric(as.character(resultsOK[,3])))
    resultsOK$geneFDR<-resultsOK$FDR<-NA
    resultsOK$FDR[!is.na(resultsOK[,"pval"])]<-p.adjust(resultsOK[!is.na(
        resultsOK[,"pval"]),"pval"], method="fdr")
    uniqGPval<-unique(resultsOK[,c("gene","genePval")])
    uniqPval<-uniqGPval[!is.na(uniqGPval[,"genePval"]),"genePval"]
    names(uniqPval)<-uniqGPval[!is.na(uniqGPval[,"genePval"]),"gene"]
    uniqPadj<-p.adjust(uniqPval, method="fdr")
    resultsOK$geneFDR<-uniqPadj[
        match(as.character(resultsOK[,"gene"]), names(uniqPadj))]
    rownames(resultsOK)<-seq_len(nrow(resultsOK))

    # Add filtered isoforms and genes mean relative expression values
    if(length(idxLowRat) > 0 ){
    geneIsoF<-isoGeneRel(object)[idxLowRat, , drop=FALSE]
    geneCountsF<-object@geneCounts[idxLowRat, , drop=FALSE]
    isoCountsF<-isoCounts(object)[idxLowRat, , drop=FALSE]
    genesF<-as.character(unique(geneIsoF[,"gene_id"]))
    resultsF<-do.call(rbind, bplapply(seq_along(genesF), function(j){
        gene<-genesF[j]
        myData<-buildData(isoCounts=isoCountsF, geneCounts=geneCountsF, 
            geneIso=geneIsoF, gene=gene, designMatrix=designMatrix,
            colName=colName)
        iso<-unique(as.character(myData[,"iso"]))
#         ratioControl<-ratioTreat<-NULL
# 
#         for(i in seq_along(iso)){
#             ratioControl<-c(ratioControl, mean(myData[myData[,colName]== 
#                 contrast[1] & myData[, "iso"] == iso[i], "counts"]/ myData[
#                 myData[,colName]==contrast[1] & myData[, "iso"] == iso[i],
#                 "all"]))
#             ratioTreat<-c(ratioTreat, mean(myData[myData[,colName]==
#                 contrast[2]& myData[, "iso"] == iso[i], "counts"]/myData[
#                 myData[,colName]==contrast[2]& myData[, "iso"] == iso[i], 
#                 "all"] ))
#         }
#
        ratios<-as.data.frame(do.call(rbind, lapply(seq_along(iso), 
            function(i){
            ratioControl<-mean(myData[myData[,colName]== contrast[1] & myData[,
                "iso"] == iso[i], "counts"]/ myData[myData[,
                colName]==contrast[1] & myData[, "iso"] == iso[i],"all"],
                na.rm=TRUE)
            ratioTreat<-mean(myData[myData[,colName]==contrast[2]& myData[,
                "iso"] == iso[i], "counts"]/myData[myData[,
                colName]==contrast[2]& myData[, "iso"] == iso[i], "all"], 
                na.rm=TRUE)
            return(c(ratioControl=ratioControl, ratioTreat=ratioTreat))
            
        })))
# testW<-data.frame(ratioControl=ratiosratioControl, ratioTreat=ratioTreat,
#     odd=NA, stat=NA, pval=NA) 
# testW<-as.data.frame(ratios)
# testW$odd<-testW$stat<-test$pval<-NA
        isoRes<-cbind(iso=iso,gene=gene, ratioControl=ratios[,"ratioControl"], 
            ratioTreat=ratios[,"ratioTreat"])
        colnames(isoRes)[c(3,4)]<-paste("ratio", contrast, sep="_")
        return(isoRes)}, BPPARAM=BPPARAM))
    resultsF<-as.data.frame(resultsF)
    for(i in c(3,4)){
        resultsF[,i]<-as.numeric(as.character(resultsF[,i]))
    }
    resultsF$theta<-NA
    resultsF$stat<-resultsF$odd<-NA
    resultsF$FDR<-resultsF$geneFDR<-resultsF$genePval<-resultsF$pval<-NA
    #join
    for(i in c(1,2)){
        resultsOK[,i]<-as.character(resultsOK[,i])            
        resultsF[,i]<-as.character(resultsF[,i])
    }
    resultsOK[,3:ncol(resultsOK)]<-vapply(3:ncol(resultsOK), function(i){
        return(as.numeric(as.character(resultsOK[,i])))
    }, as.numeric(as.character(resultsOK[,3])))
    resultsAll<-rbind(resultsOK, resultsF)
    }else{resultsAll<-resultsOK}
    rownames(resultsAll)<-as.character(resultsAll[,"iso"])

    iso_cm<-isoCounts(object)
    resultsAll<-resultsAll[rownames(iso_cm),]
    rownames(resultsAll)<-seq_len(nrow(resultsAll))
    
    resultsAll$genePval<-uniqPval[match(as.character(resultsAll[,"gene"]),
        names(uniqPval))]
    
    resultsAll$geneFDR<-uniqPadj[match(as.character(resultsAll[,"gene"]),
                                        names(uniqPadj))]
    
    genes<-unique(resultsAll[,"gene"])
    disp<-do.call(c,lapply(genes, function(gene){
        if(any(!is.na(resultsAll[resultsAll[,"gene"]==gene,"theta"]))){
            disp<-unique(resultsAll[resultsAll[,"gene"]==gene & !is.na(
                resultsAll[, "theta"]),"theta"])
        }else{
            disp<-NA}
        return(disp)
    }))
    names(disp)<-genes
    resultsAll<-resultsAll[, colnames(resultsAll) != "theta"]
    NBSpliceResults<-NBSpliceRes(resultsAll, idxLowRat, contrast, disp)
    return(NBSpliceResults)       
})
