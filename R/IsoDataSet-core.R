#'@title
#'GLM fitting and hypothesis testing
#'@description
#'\code{buildData} is an R function to build the data
#'structure useful to fit a negative binomial model at the gene level.
#'\code{fitModel} is an R function to fit a negative binomial generalized
#'linear model in order to evaluate changes in the isoform ratio. Each isoform
#'is considered as a factor. The model incorporates one factor related to one
#'experimental condition with two levels.
#'@param isoCounts Matrix having the expression counts at the isoform level. 
#'Isoforms must be in rows and samples in columns. Rownames and colnames must 
#'be defined with isoform and samples names, respectively.
#'@param geneCounts Matrix having the expression counts at the gene level. 
#'Genes must be in rows and samples in columns. Rownames and colnames must 
#'be defined with gene and sample names, respectively.
#'@param geneIso Data.frame containing the relationship between isoforms and 
#'genes. It must contain two columns, named as 'gene_id' and 'isoform_id'.
#'Its isoforms should be the same specified in the isoCounts matrix.
#'@param gene Character indicating the name of the gene to be analyzed
#'@param designMatrix Data.frame specifying metadata related to the
#'experiment. Its rows must be the samples and experimental factors should be
#'arranged on its columns.
#'@param colName Character indicating the name of the column in the design 
#'matrix to be considered for mean expression calculations per experimental
#'condition and differential expression test.
#'@return A data.frame ready to use by the model fitting function
#'@include IsoDataSet-buildLowExpIdx.R
#'@docType methods
#'@name buildData
#'@rdname IsoDataSet-core
#'@aliases buildData-methods
#'@seealso \code{\link{IsoDataSet}}
#'@note see full example in \code{\link{IsoDataSet-class}}
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#'Fernandez \email{efernandez@@bdmg.com.ar}
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
buildData<-function(isoCounts, geneCounts, geneIso, gene, designMatrix,
    colName){
    iso_cm<-cbind(isoCounts, geneCounts)
    if(!any(colnames(designMatrix)==colName)){
        stop(paste("The design matrix should contain a column called", colName,
            sep=" "))
    }
    iso_cm_gen<-iso_cm[which(geneIso[,"gene_id"] == gene),,drop=FALSE]
    iso<-rownames(iso_cm_gen)
    data<-data.frame(samples=rep(rownames(designMatrix), length(iso)), 
        condition=rep(designMatrix[,colName], length(iso)))
    for (i in seq_along(iso)){
        data[((i+ (i-1)*(nrow(designMatrix)-1)):(i+ (i-1)*(nrow(
            designMatrix)-1)+nrow(designMatrix)-1)), "iso"]<-iso[i]
    }
    data[,"counts"]<-as.numeric(t(iso_cm_gen[,rownames(designMatrix)]))
    data[,"all"]<-as.numeric(t(iso_cm_gen[,paste(rownames(designMatrix),
        "_All", sep="")]))
    data[,"iso"]<-factor(as.character(data[,"iso"]), levels=unique( 
        data[,"iso"]))
    colnames(data)[colnames(data)=="condition"]<-colName
    return(data)
}
#'@param myData Data.frame containing the expression matrix at isoform 
#'levels, it means, isoform in rows and samples in columns. It is obtained 
#'using the buildData NBSplice method.
#'@param formula Object with the formula of the GLM.
#'@param test Character indicating the name of the distribution to assume for
#'linear hypothesis statistic. Could be "F" or "chisq".
#'@param contrast Character vector with the names of the two levels of the
#'experimental factor to be contrasted.
#'@return A data.frame summarizing gene results.
#'@docType methods
#'@name fitModel
#'@rdname IsoDataSet-core
#'@inheritParams buildData
#'@import MASS
#'@importFrom MASS glm.nb
#'@import stats
#'@importFrom stats drop1
#'@importFrom stats coefficients
#'@importFrom stats glm.control
#'@importFrom stats residuals
#'@import car
#'@importFrom car lht
#'@import mppa
#'@importFrom mppa simes.test
#'@aliases fitModel-methods
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#' Fernandez \email{efernandez@bdmg.com.ar}
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
fitModel<-function(myData, gene, formula, colName, test=c("F", "Chisq"), 
    contrast){
    iso<-levels(myData[, "iso"])
    myData[,"counts"]<-round(myData[,"counts"])
    myData[,"all"]<-round(myData[,"all"])
    myData<-myData[as.character(myData[, "condition"]) %in% contrast, ]
    myData[,"condition"]<-factor(as.character(myData[, "condition"]), 
        levels=contrast)
    if(length(levels(myData[,"iso"]))> 1 & all(myData[,"all"] > 0) & !(all(
        myData[,"counts"] ==(myData[,"all"] -myData[,"counts"])))){
#             modl<-tryCatch(glm.nb(formula,
            modl<-suppressWarnings(
                tryCatch(glm.nb(counts~condition+iso+condition:iso, 
                offset=log(myData[,"all"]), 
                link="log", data=myData, control=glm.control()),
                error=function(cond){FALSE}))
            if(is.logical(modl)){   
                converged<-modl 
            }else{
                converged<-modl$converged
            }
    }else{
        if(all(myData[,"all"] > 0) & sum(myData[,"all"]) > 0 & !(all(
            myData[,"counts"] == (myData[,"all"] -myData[,"counts"])))){
                form<-as.formula(paste("counts~", colName))
#                 modl<-tryCatch( glm.nb(form, offset=log(
                modl<-suppressWarnings(
                    tryCatch( glm.nb(counts~condition, offset=log(

                    myData[,"all"]), link="log", data=myData, 
                    control=glm.control()), error=function(cond){FALSE}))
            if(is.logical(modl)){
                converged<-modl
            }else{
                converged<-modl$converged
            }
        }else{
            converged<-FALSE
        }
    }
    if(converged){
        if(modl$df.residual  >1 & sum(residuals(modl, "pearson")^2)>0){
#             Ftest<-drop1(modl,test=test)
#             FBeta<-Ftest[,4][2]
            beta<-coefficients(modl)
            sigma2<-sum(residuals(modl, "pearson")^2)/modl$df.residual
            x<-model.matrix(modl)
            theta<-modl$theta
            testW<-do.call(rbind,lapply(seq_along(iso), function(k){
                isof<-iso[k]
                coef<-which(names(beta) %in% c(paste(colName, contrast[2],
                    sep=""), paste(colName, contrast[2], ":iso", isof, 
                    sep="")))
                if(length(coef)== 0){
                    coef<-which(names(beta) == paste(colName, contrast[2],
                        sep=""))
                }
                H<-matrix(0,nrow=1, ncol=length(beta))
                H[1,coef]<-1
                aux<-tryCatch({lht(modl, H, test=test)[2,] }, error= function(
                    cond){
                        aux<-data.frame(stat=NA, prob=NA)
                        colnames(aux)<-c("stat", "prob")
                        return(aux)
                })
                if(any(colnames(aux) %in% c("F", "Chisq")) | any(colnames(aux)
                    %in% c("Pr(>F)","Pr(>Chisq)" ))){
                        colnames(aux)[colnames(aux) %in% c("F", 
                            "Chisq")]<-"stat"
                        colnames(aux)[colnames(aux) %in% c("Pr(>F)",
                            "Pr(>Chisq)")]<-"prob"
                }
                ratioControl<-exp(sum(as.numeric(beta[c(1,which(names(
                    beta) %in% paste("iso", isof, sep="")))])))
                ratioTreat<-exp(sum(c(as.numeric(beta[c(1,which(names(
                    beta) %in% paste("iso", isof, sep="")))]), beta[coef])))
                return(c(ratioControl=ratioControl, ratioTreat=ratioTreat, 
                    odd=sum(beta[coef]), stat=aux[,"stat"], pval=aux[,"prob"]))
            }))
            genePval<-simes.test(as.numeric(testW[,"pval"]))
        }else{
            testW<-data.frame(ratioControl=NA, ratioTreat=NA, odd=NA, stat=NA, 
                pval=NA)
#             theta<-sigma2<-FBeta<-genePval<-NA
            theta<-genePval<-NA
        }
    }else{
        ratioControl<-ratioTreat<-NULL
        for(i in seq_along(iso)){
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
    }
    isoRes<-cbind(iso=iso,gene=gene, ratioControl=testW[,"ratioControl"], 
        ratioTreat=testW[,"ratioTreat"],theta=theta, odd=testW[,"odd"], 
        stat=testW[,"stat"], pval=testW[,"pval"], genePval=genePval)

    colnames(isoRes)[3:4]<-paste("ratio", contrast, sep="_")
    return(isoRes)
}
