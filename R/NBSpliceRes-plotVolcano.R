#'@title
#'Method to obtain the Volcano Plot. 
#'@description
#'\code{plotVolcano} returns a generalized volcano plot where the x-axis
#'represents the difference between isoform's relative expression in the 
#'contrasted conditions. Isoforms are colored according to the differential
#'splicing status of the gene where they come from.
#'@param myNBRes NBSpliceRes class object.
#'@param adjusted Logical indicating if adjusted p values should be used.
#'@param p.value Numeric value between 0 and 1 giving the required family-wise
#'error rate or false discovery rate.
#'
#'@return A ggplot object.
#'
#'@include NBSpliceRes-plotRatiosDisp.R
#'@exportMethod plotVolcano
#'@docType methods
#'@name plotVolcano
#'@rdname NBSpliceRes-plotVolcano
#'@import ggplot2
#'@import reshape2
#'@importFrom reshape2 melt
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_point
#'@aliases plotVolcano-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'g<-plotVolcano(myDSResults)
#'if(interactive()){
#'g
#'}
setGeneric(name="plotVolcano", def=function(myNBRes, adjusted=TRUE, 
    p.value=0.05){
        standardGeneric("plotVolcano")
})
#'@name plotVolcano
#'@rdname NBSpliceRes-plotVolcano
#'@aliases plotVolcano,NBSpliceRes-method
#'@inheritParams plotVolcano
setMethod(f="plotVolcano", signature="NBSpliceRes", definition=function(
    myNBRes, adjusted=TRUE, p.value=0.05){
        DSDF<-results(myNBRes)
        DSDF<-DSDF[!is.na(DSDF[, "pval"]),]
        condVars<-colnames(DSDF)[grep("ratio", colnames(DSDF))]
        colnames(DSDF)[grep("ratio", colnames(DSDF))]<-c("x", "y")
        DSDF$ratDif<-DSDF[,"y"]-DSDF[,"x"]
        ratDif<-FDR<-DiffSpl<-pval<-NULL
        if(adjusted){
            DSDF[,"DiffSpl"]<-DSDF[,"geneFDR"] < p.value & !is.na( DSDF[,
                "geneFDR"])
            g<-ggplot(DSDF, aes(x=ratDif, y=-log10(FDR), color=DiffSpl))+
            geom_point()+labs(x=paste(condVars[2], "-", condVars[1], sep=""))
        }else{ 
            DSDF[,"DiffSpl"]<-DSDF[,"genePval"] < p.value & !is.na( DSDF[,
                "genePval"])
            g<-ggplot(DSDF, aes(x=ratDif, y=-log10(pval), color=DiffSpl))+
            geom_point()+labs(x=paste(condVars[2], "-", condVars[1], sep=""))
        }
        g<-g+theme(panel.background=element_rect(fill="white", color="black"),
        legend.key=element_rect( fill="white", color="white"), 
        panel.grid.major=element_line(color="lightgoldenrod3"),
        panel.grid.minor=element_line(color="tomato", linetype = "dashed"),
        legend.text = element_text(size = 12),legend.title = element_text(
        size = 14, face = "bold"), axis.text=element_text(size=12), 
        plot.title=element_text(size=14, hjust=0.5), axis.title=element_text(
        size=14))+scale_color_manual(values=c("TRUE"="green", "FALSE"="red"))
        return(g)
})