#'@title
#'Plot the dispersion of isoform's relative expression for significant genes. 
#'@description
#'\code{plotRatiosDisp} returns a ggplot object with the scatter plot of 
#'isoform's relative expression in the contrasted conditions. Isoforms are 
#'colored according to the differential splicing status of the gene where they
#'come from.
#'@param myNBRes NBSpliceRes class object.
#'@param adjusted Logical indicating if adjusted p values should be used.
#'@param p.value Numeric value between 0 and 1 giving the required family-wise
#'error rate or false discovery rate.
#'
#'@return A ggplot object.
#'
#'@include NBSpliceRes-getGeneResults.R
#'@exportMethod plotRatiosDisp
#'@docType methods
#'@name plotRatiosDisp
#'@rdname NBSpliceRes-plotRatiosDisp
#'@import ggplot2
#'@import reshape2
#'@importFrom reshape2 melt
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_point
#'@aliases plotRatiosDisp-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'
#'g<-plotRatiosDisp(myDSResults)
#'if(interactive()){
#'g
#'}
setGeneric(name="plotRatiosDisp", def=function(myNBRes, adjusted=TRUE, 
    p.value=0.05){
        standardGeneric("plotRatiosDisp")
})
#'@name plotRatiosDisp
#'@rdname NBSpliceRes-plotRatiosDisp
#'@aliases plotRatiosDisp,NBSpliceRes-method
#'@inheritParams plotRatiosDisp
setMethod(f="plotRatiosDisp", signature="NBSpliceRes", definition=function(
    myNBRes, adjusted=TRUE, p.value=0.05){
        DSDF<-results(myNBRes)
        DSDF<-DSDF[!is.na(DSDF[, "pval"]),]
        condVars<-colnames(DSDF)[grep("ratio", colnames(DSDF))]
        colnames(DSDF)[grep("ratio", colnames(DSDF))]<-c("x", "y")
        
        if(adjusted){
            DSDF[,"DiffSpl"]<-DSDF[,"geneFDR"] < p.value & !is.na( DSDF[,
                "geneFDR"])
        }else{        
            DSDF[,"DiffSpl"]<-DSDF[,"genePval"] < p.value & !is.na( DSDF[,
                "genePval"])
        }
        x<-y<-DiffSpl<-NULL
        g<-ggplot(DSDF, aes(x=x, y=y, color=DiffSpl))+ geom_point()+labs(x=
            condVars[1], y=condVars[2])

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