library(RUnit);
library(NBSplice);

rOpts <- getOption("RUnit");
rOpts$silent <- !FALSE;
options("RUnit"=rOpts);

require(BiocParallel);
require(ggplot2)
if (.Platform$OS.type == "unix") {
    bp_param <- MulticoreParam(workers=1);
} else if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers=1);
}
register(bp_param);

##-----------------------------------------------------------------------------
##Constructor-class Tests 
##-----------------------------------------------------------------------------

test_IsoDataSet<-function(){
    checkTrue(validObject(IsoDataSet()), 
              msg="IsoDataSet contstructor works without any parameter: OK.") 
}

test_NBSpliceRes<-function(){
    checkTrue(validObject(NBSpliceRes()), 
              msg="NBSpliceRes contstructor works without any parameter: OK.") 
}

##Build IsoDataSet object with user data
test_IsoDataSetWithData<-function(){
    data(isoCounts, package="NBSplice")
    data(designMatrix, package="NBSplice")
    data(geneIso, package="NBSplice")
    colName<-"condition"    
    checkEquals(ncol(isoCounts), nrow(designMatrix),
                msg="isoCounts and designMatrix consistency: OK.")
    checkEquals(nrow(isoCounts), nrow(geneIso),
                msg="isoCounts and geneIso consistency: OK.")
    checkTrue(all(rownames(isoCounts)==rownames(geneIso)),
                msg="isoCounts and geneIso order: OK.")
    checkTrue(any(colnames(designMatrix)==colName),
                msg="designMatrix and colName consistency: OK.")
    
    object<-IsoDataSet(isoCounts = isoCounts, experimentData = designMatrix,
                colName = colName, geneIso = geneIso)
    # Checking slots
    #counts
    checkTrue(is.matrix(isoCounts(object)), 
                msg="IsoDataSet isoCounts slot type: OK.")
    #geneCounts
    checkTrue(is.matrix(geneCounts(object)),  
                msg="IsoDataSet geneCounts slot type: OK.")
    #colData and design
    checkTrue(is.data.frame(expData(object)),  
                msg="IsoDataSet colData slot type: OK.")
    checkEquals(class(designFormula(object)),"formula",  
                msg="IsoDataSet design slot type: OK.")
    #lowExpIndex
    checkTrue(is.integer(lowExpIndex(object)), 
                msg="LowExpIndex slot type: OK.")
}

##-----------------------------------------------------------------------------
##Getters Tests 
##-----------------------------------------------------------------------------
test_isoCounts<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkTrue(is.matrix(isoCounts(myIsoDataSet)), 
            msg="IsoDataSet isoCounts slot type: OK.")
}
test_geneCounts<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkTrue(is.matrix(geneCounts(myIsoDataSet)),  
            msg="IsoDataSet geneCounts slot type: OK.")
}
test_expData<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkTrue(is.data.frame(expData(myIsoDataSet)),  
            msg="IsoDataSet colData slot type: OK.")
}
test_designFormula<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkEquals(class(designFormula(myIsoDataSet)),"formula",  
            msg="IsoDataSet design slot type: OK.")
}
test_isoGeneRel<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkTrue(is.data.frame(isoGeneRel(myIsoDataSet)),  
            msg="IsoDataSet isoGeneRel slot type: OK.")
}
test_lowExpIndex<-function(){
    data(myIsoDataSet, package="NBSplice")
    checkTrue(is.integer(lowExpIndex(myIsoDataSet)), 
            msg="LowExpIndex slot type: OK.")
}
test_results<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.data.frame(results(myDSResults)),
            msg="Results slot type: OK.")
    checkTrue(nrow(results(myDSResults))>0 & ncol(results(myDSResults))==10 ,
            msg="Results slot dimention: OK.")
    checkTrue(all(c("iso", "gene", "FDR", "geneFDR") %in% colnames(results(
        myDSResults))), msg="Results slot columns: OK.")
    checkEquals(dim(results(NBSpliceRes())),c(0,0),
            msg = "Results slot from an zero length object: OK.")
}
test_contrast<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.character(contrast(myDSResults)),
            msg="Contrast slot type: OK.")
    checkTrue(all(contrast(myDSResults) %in% c("Normal", "Tumor")),
            msg="Contrast slot: OK.")
    checkEquals(length(contrast(NBSpliceRes())),0,
            msg = "Contrast slot from an zero length object: OK.")
}
test_disp<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.numeric(disp(myDSResults)),
            msg="Dispersion slot type: OK.")
    checkEquals(length(disp(NBSpliceRes())),0,
            msg = "Dispersion slot from an zero length object: OK.")
}
test_GetDSGenes<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.character(GetDSGenes(myDSResults)), 
            msg="Getting DS genes type: OK.")
    checkEquals(length(GetDSGenes(NBSpliceRes())),0,
            msg = "Getting DS genes from an zero length object: OK.")
}
test_GetDSResults<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.data.frame(GetDSResults(myDSResults)), 
            msg="Getting DS results type: OK.")
    checkEquals(dim(GetDSResults(myDSResults)),c(110,10),
            msg = "Getting DS results dimension: OK.")
    checkEquals(dim(GetDSResults(NBSpliceRes())),c(0,0),
            msg = "Getting DS results dimension from an zero length object: 
                OK.")
}
test_GetGeneResults<-function(){
    data(myDSResults, package="NBSplice")
    checkTrue(is.data.frame(GetGeneResults(myDSResults, 
        gene = "ENSG00000005889")), msg="Getting DS gene results type: OK.")
    checkEquals(dim(GetGeneResults(myDSResults, gene = "ENSG00000005889")),
                c(3,10), msg = "Getting DS results dimension: OK.")
}

##-----------------------------------------------------------------------------
##NBSplice core function Tests 
##-----------------------------------------------------------------------------
test_buildLowExpIdx<-function(){
    data(myIsoDataSet, package="NBSplice")
    aux<-buildLowExpIdx(myIsoDataSet)
    checkTrue(length(lowExpIndex(aux)) > length(lowExpIndex(myIsoDataSet)), 
            msg="buildLowExp method performance: OK.")
    checkEqualsNumeric(round(length(lowExpIndex(aux))), 2224, 
            msg="buildLowExp length: OK.")
    checkTrue(is.integer(lowExpIndex(aux)), msg="buildLowExp object type: OK.")
    aux2<-buildLowExpIdx(myIsoDataSet, ratioThres = 1)
    checkEquals(length(lowExpIndex(aux2)), 3144, 
            msg="buildLowExp method performance extreme filtering value: OK.")
    
}
test_NBTest<-function(){
    data(geneIso, package="NBSplice")
    data(isoCounts, package="NBSplice")
    data(designMatrix, package="NBSplice")
    colName<-"condition"
    test<-"F"
    GIR<-geneIso[geneIso[,"gene_id"]=="ENSG00000002016",]
    IC<-isoCounts[GIR[,"isoform_id"],]
    obj<-IsoDataSet(isoCounts = IC, experimentData=designMatrix, colName, 
        geneIso=GIR)
    obj2<-buildLowExpIdx(obj)
    res<-suppressWarnings(NBTest(object=obj2, colName, test))
    checkEquals(class(res)[1], "NBSpliceRes", 
            msg="Checking NBTest results type: OK.")
    checkEquals(dim(results(res, filter=FALSE)), c(20,10), 
            msg="Checking NBTest results dimention: OK.")
    checkEquals(dim(results(res, filter=TRUE)), c(1,10), 
            msg="Checking NBTest results dimention 2: OK.")
    checkTrue(all(!is.na(results(res, filter=TRUE)[, "geneFDR"])),  
            msg="Checking NBTest results computation: OK.")
}
##-----------------------------------------------------------------------------
##NBSplice plotting methods tests 
##-----------------------------------------------------------------------------
test_plotGeneResults<-function(){
    data(myDSResults, package="NBSplice")
    g<-plotGeneResults(myDSResults, gene="ENSG00000002016")
    checkTrue(is.ggplot(g), msg = "Checking plotGeneResults class: OK.")
}
test_plotRatiosDisp<-function(){
    data(myDSResults, package="NBSplice")
    g<-plotRatiosDisp(myDSResults)
    checkTrue(is.ggplot(g), msg = "Checking plotRatiosDisp class: OK.")
}
test_plotVolcano<-function(){
    data(myDSResults, package="NBSplice")
    g<-plotVolcano(myDSResults)
    checkTrue(is.ggplot(g),msg = "Checking plotVolcano class: OK.")
}

