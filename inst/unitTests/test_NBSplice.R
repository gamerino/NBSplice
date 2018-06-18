library(RUnit);
library(NBSplice);

rOpts <- getOption("RUnit");
rOpts$silent <- !FALSE;
options("RUnit"=rOpts);

require(BiocParallel);
if (.Platform$OS.type == "unix") {
    bp_param <- MulticoreParam(workers=1);
} else if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers=1);
}
register(bp_param);

data(myIsoDataSet, package="NBSplice")
data(myNBSpliceRes, package="NBSplice")
data(isoCounts, package="NBSplice")
data(designMatrix, package="NBSplice")
data(geneIso, package="NBSplice")
colName<-"condition"    
test<-"F"
##-----------------------------------------------------------------------------
##Constructor-class Tests 
##-----------------------------------------------------------------------------

test_IsoDataSet(){
    checkTrue(validObject(IsoDataSet()), 
              msg="IsoDataSet contstructor works without any parameter: OK.") 
}

test_NBSpliceRes(){
    checkTrue(validObject(NBSpliceRes()), 
              msg="NBSpliceRes contstructor works without any parameter: OK.") 
}

##Build IsoDataSet object with user data
test_IsoDataSetWithData<-function(){
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
    checkTrue(is.matrix(isoCounts(object)), 
            msg="IsoDataSet isoCounts slot type: OK.")
}
test_geneCounts<-function(){
    checkTrue(is.matrix(geneCounts(object)),  
            msg="IsoDataSet geneCounts slot type: OK.")
}
test_expData<-function(){
    checkTrue(is.data.frame(expData(object)),  
            msg="IsoDataSet colData slot type: OK.")
}
test_designFormula<-function(){
    # loading data matrices
    checkEquals(class(designFormula(object)),"formula",  
            msg="IsoDataSet design slot type: OK.")
}
test_isoGeneRel<-function(){
    # loading data matrices
    checkTrue(is.data.frame(isoGeneRel(object)),  
            msg="IsoDataSet isoGeneRel slot type: OK.")
}
test_lowExpIndex<-function(){
    checkTrue(is.integer(lowExpIndex(object)), 
            msg="LowExpIndex slot type: OK.")
}
test_results<-function(){
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
    checkTrue(is.character(contrast(myDSResults)),
            msg="Contrast slot type: OK.")
    checkTrue(all(contrast(myDSResults) %in% c("Normal", "Tumor")),
            msg="Contrast slot: OK.")
    checkEquals(length(contrast(NBSpliceRes())),0,
            msg = "Contrast slot from an zero length object: OK.")
}
test_disp<-function(){
    checkTrue(is.numeric(disp(myDSResults)),
            msg="Dispersion slot type: OK.")
    checkEquals(length(disp(NBSpliceRes())),0,
            msg = "Dispersion slot from an zero length object: OK.")
}
test_GetDSGenes<-function(){
    checkTrue(is.character(GetDSGenes(myDSResults)), 
            msg="Getting DS genes type: OK.")
    checkEquals(length(GetDSGenes(NBSpliceRes())),0,
            msg = "Getting DS genes from an zero length object: OK.")
}
test_GetDSResults<-function(){
    checkTrue(is.data.frame(GetDSResults(myDSResults)), 
            msg="Getting DS results type: OK.")
    checkEquals(dim(GetDSResults(myDSResults)),c(126,10),
            msg = "Getting DS results dimension: OK.")
    checkEquals(dim(GetDSResults(NBSpliceRes())),c(0,0),
            msg = "Getting DS results dimension from an zero length object: 
                OK.")
}
test_GetGeneResults<-function(){
    checkTrue(is.data.frame(GetGeneResults(myDSResults, 
        gene = "ENSG00000010256")), msg="Getting DS gene results type: OK.")
    checkEquals(dim(GetGeneResults(myDSResults, gene = "ENSG00000010256")),
                c(9,10), msg = "Getting DS results dimension: OK.")
}

##-----------------------------------------------------------------------------
##NBSplice core function Tests 
##-----------------------------------------------------------------------------
test_buildLowExpIdx<-function(){
    aux<-buildLowExpIdx(myIsoDataSet)
    checkTrue(lowExpIndex(aux)!=lowExpIndex(myIsoDataSet), 
            msg="buildLowExp method performance: OK.")
    checkEquals(length(lowExpIndex(aux)), 2146, 
            msg="buildLowExp length: OK.")
    checkTrue(is.integer(lowExpIndex(aux)), msg="buildLowExp object type: OK.")
    aux2<-buildLowExpIdx(myIsoDataSet, ratioThres = 1)
    checkEquals(length(lowExpIndex(aux2)), 3024, 
            msg="buildLowExp method performance extreme filtering value: OK.")
    
}
test_NBTest<-function(){
    GIR<-geneIso[geneIso[,"gene_id"]=="ENSG00000006704",]
    IC<-isoCounts[GIR[,"isoform_id"],]
    obj<-IsoDataSet(isoCounts = IC, experimentData=designMatrix, colName, 
        geneIso=GIR)
    obj2<-buildLowExpIdx(obj)
    res<-suppressWarnings(NBTest(obj2, colName, test))
    checkEquals(class(res)[1], "NBSpliceRes", 
            msg="Checking NBTest results type: OK.")
    checkTrue(dim(results(res, filter=FALSE)), c(3,10), 
            msg="Checking NBTest results dimention: OK.")
    checkTrue(dim(results(res, filter=TRUE)), c(0,10), 
            msg="Checking NBTest results dimention 2: OK.")
    checkTrue(all(!is.na(results(res, filter=F)[, "geneFDR"])),  
            msg="Checking NBTest results computation: OK.")
}
##-----------------------------------------------------------------------------
##NBSplice plotting methods tests 
##-----------------------------------------------------------------------------
test_plotGeneResults<-function(){
    g<-plotGeneResults(myDSResults, gene="ENSG00000006704")
    checkTrue(is.ggplot(g), msg = "Checking plotGeneResults class: OK.")
}
test_plotRatiosDisp<-function(){
    g<-plotRatiosDisp(myDSResults)
    checkTrue(is.ggplot(g), msg = "Checking plotRatiosDisp class: OK.")
}
test_plotVolcano<-function(){
    g<-plotVolcano(myDSResults)
    checkTrue(is.ggplot(g),msg = "Checking plotVolcano class: OK.")
}
