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
    test=c("F", "Chisq"), contrast=c(levels(expData(object)[,colName])[1:2]), 
    BPPARAM=bpparam()){
    standardGeneric("NBTest")
})
#'@name NBTest
#'@rdname IsoDataSet-NBTest
#'@aliases NBTest,IsoDataSet-method
#'@inheritParams NBTest
setMethod(f="NBTest", signature=signature(object="IsoDataSet"),
    definition=function(object, colName="condition", test=c("F", "Chisq"), 
    contrast=c(levels(expData(object)[,colName])[1:2]), BPPARAM=bpparam()){
    idxLowRat<-lowExpIndex(object)
    if(length(idxLowRat)==0){
        object<-buildLowExpIdx(object)
        idxLowRat<-lowExpIndex(object)
    }
    designMatrix<-expData(object)
    geneIso<-isoGeneRel(object)[-idxLowRat,]
    geneCounts<-object@geneCounts[-idxLowRat,]
    isoCounts<-isoCounts(object)[-idxLowRat,]
    genes<-as.character(unique(geneIso[,"gene_id"]))
    
    resultsOK<-do.call(rbind, bplapply(1:length(genes), function(j){
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
    for(i in 3:ncol(resultsOK)){
        resultsOK[,i]<-as.numeric(as.character(resultsOK[,i]))
    }
    resultsOK$geneFDR<-resultsOK$FDR<-NA
    resultsOK$FDR[!is.na(resultsOK[,"pval"])]<-p.adjust(resultsOK[!is.na(
        resultsOK[,"pval"]),"pval"], method="fdr")
    resultsOK$geneFDR[!is.na(resultsOK[,"genePval"])]<-p.adjust(resultsOK[
        !is.na(resultsOK[,"genePval"]),"genePval"], method="fdr")
    rownames(resultsOK)<-1:nrow(resultsOK)
    
    # Add filtered isoforma and genes mean relative expression values
    geneIsoF<-isoGeneRel(object)[idxLowRat,]
    geneCountsF<-object@geneCounts[idxLowRat,]
    isoCountsF<-isoCounts(object)[idxLowRat,]
    genesF<-as.character(unique(geneIsoF[,"gene_id"]))
    resultsF<-do.call(rbind, bplapply(1:length(genesF), function(j){
        gene<-genesF[j]
        myData<-buildData(isoCounts=isoCountsF, geneCounts=geneCountsF, 
            geneIso=geneIsoF, gene=gene, designMatrix=designMatrix,
            colName=colName)
        iso<-unique(as.character(myData[,"iso"]))
        ratioControl<-ratioTreat<-NULL
        
        for(i in 1:length(iso)){
            ratioControl<-c(ratioControl, mean(myData[myData[,colName]== 
                contrast[1] & myData[, "iso"] == iso[i], "counts"]/ myData[
                myData[,colName]==contrast[1] & myData[, "iso"] == iso[i],
                "all"]))
            ratioTreat<-c(ratioTreat, mean(myData[myData[,colName]==
                contrast[2]& myData[, "iso"] == iso[i], "counts"]/myData[
                myData[,colName]==contrast[2]& myData[, "iso"] == iso[i], 
                "all"] ))
        }
        testW<-data.frame(ratioControl=ratioControl, ratioTreat=ratioTreat,
            odd=NA, stat=NA, pval=NA) 
#         theta<-sigma2<-FBeta<-genePval<-NA
        theta<-genePval<-NA
        isoRes<-cbind(iso=iso,gene=gene, ratioControl=testW[,"ratioControl"], 
            ratioTreat=testW[,"ratioTreat"])
        colnames(isoRes)[3:4]<-paste("ratio", contrast, sep="_")
        return(isoRes)}, BPPARAM=BPPARAM))
    resultsF<-as.data.frame(resultsF)
    for(i in 3:4){
        resultsF[,i]<-as.numeric(as.character(resultsF[,i]))
    }
    resultsF$theta<-NA
    resultsF$stat<-resultsF$odd<-NA
    resultsF$FDR<-resultsF$geneFDR<-resultsF$genePval<-resultsF$pval<-NA
    #join
    for(i in 1:2){
        resultsOK[,i]<-as.character(resultsOK[,i])            
        resultsF[,i]<-as.character(resultsF[,i])
    }
    resultsAll<-rbind(resultsOK, resultsF)
    
    rownames(resultsAll)<-as.character(resultsAll[,"iso"])
    
    iso_cm<-isoCounts(object)
    resultsAll<-resultsAll[rownames(iso_cm),]
    rownames(resultsAll)<-1:nrow(resultsAll)
    
    genes<-unique(resultsAll[,"gene"])
    disp<-do.call(c,lapply(genes, function(gene){
        if(any(!is.na(resultsAll[resultsAll[,"gene"]==gene,"theta"]))){
            disp<-unique(resultsAll[resultsAll[,"gene"]==gene & !is.na(resultsAll[,
                "theta"]),"theta"])
        }else{
            disp<-NA}
        return(disp)
    }))
    names(disp)<-genes
    resultsAll<-resultsAll[, colnames(resultsAll) != "theta"]
    NBSpliceResults<-NBSpliceRes(resultsAll, idxLowRat, contrast, disp)
    return(NBSpliceResults)       
})