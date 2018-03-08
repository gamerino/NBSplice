## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE, 
    error=FALSE, warning=TRUE)

## ----installDev, eval=FALSE----------------------------------------------
## install.packages("devtools")

## ----libDev, eval=FALSE--------------------------------------------------
## library(devtools)

## ----install, eval=FALSE-------------------------------------------------
## install_github("gamerino/NBSplice")

## ----echo=FALSE----------------------------------------------------------
library(NBSplice)

## ----expMat, eval=TRUE---------------------------------------------------
data(isoCounts, package="NBSplice")
head(isoCounts)
dim(isoCounts)

## ----geneIsoMat, eval=TRUE-----------------------------------------------
data(geneIso, package="NBSplice")
head(geneIso)
dim(geneIso)

## ----designMat, eval=TRUE------------------------------------------------
data(designMatrix, package="NBSplice")
head(designMatrix)
dim(designMatrix)

## ----colData, eval=TRUE--------------------------------------------------
colName<-"condition"
levels(designMatrix[,colName])

## ----object loading, echo=FALSE, eval=TRUE, results="hide"---------------
data(myIsoDataSet, package="NBSplice")

## ----objectBuild, eval=FALSE---------------------------------------------
## myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)

## ----objectInsp, eval=TRUE-----------------------------------------------
show(myIsoDataSet)

head(isoCounts(myIsoDataSet))

## ----lowExp, eval=TRUE---------------------------------------------------
myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01, countThres = 1)

## ----NBTest, eval=FALSE--------------------------------------------------
## myDSResults<-NBTest(myIsoDataSet, colName, test="F")

## ----object loading2, echo=FALSE, eval=TRUE, results="hide"--------------
data(myDSResults, package="NBSplice")

## ----NBRes, eval=TRUE----------------------------------------------------
head(lowExpIndex(myDSResults))
contrast(myDSResults)
head(results(myDSResults))

## ----lowExpRes, eval=TRUE------------------------------------------------
head(results(myDSResults, filter = FALSE))

## ----getDSGe, eval=TRUE--------------------------------------------------
mySignificantRes<-GetDSResults(myDSResults)
head(mySignificantRes)
dim(mySignificantRes)
myDSGenes<-GetDSGenes(myDSResults)
head(myDSGenes)
length(myDSGenes)

## ----plotD, eval=TRUE----------------------------------------------------
plotRatiosDisp(myDSResults)

## ----plotv, eval=TRUE----------------------------------------------------
plotVolcano(myDSResults)

## ----plotGene, eval=TRUE-------------------------------------------------
## Select the first differentially spliced gene
gene<-GetDSGenes(myDSResults)[1]
GetGeneResults(myDSResults, gene)
plotGeneResults(myDSResults, gene)

## ----plotGene2, eval=TRUE------------------------------------------------
## Keeping non-reliable and non-significant isoforms
plotGeneResults(myDSResults, gene, filterLowExpIso = FALSE, filterNotSignificant = FALSE)


## ----sessionInfo---------------------------------------------------------
sessionInfo()

