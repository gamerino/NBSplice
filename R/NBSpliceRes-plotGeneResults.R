#'@title
#'Method to obtain isoform's relative expression barplot for an specific gene. 
#'@description
#'\code{plotGeneResults} returns ggplot object which illustrates the isoform's 
#'relative expression in the two contrasted conditions.
#'@param myNBRes NBSpliceRes class object.
#'@param gene Character indicating the gene name. 
#'@param filterLowExpIso Logical indicating if lower-expression isoforms should
#'be filtered out.
#'@param filterNotSignificant Logical indicating if not significant 
#'isoforms should be filtered out. 
#'@param adjusted Logical indicating if adjusted p values should be used.
#'@param p.value Numeric value between 0 and 1 giving the required family-wise
#'error rate or false discovery rate.
#'@param group Logical indicating if isoforms bars should be stacked or not
#'@return A ggplot object.
#'@include NBSpliceRes-plotVolcano.R
#'@exportMethod plotGeneResults
#'@docType methods
#'@name plotGeneResults
#'@rdname NBSpliceRes-plotGeneResults
#'@import ggplot2
#'@import reshape2
#'@importFrom reshape2 melt
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_point
#'@aliases plotGeneResults-methods
#'@seealso \code{\link{NBSpliceRes}}
#'@note see full example in \code{\link{NBSpliceRes-class}}
#'@family NBSpliceRes
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A. 
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@examples
#'data(myDSResults, package="NBSplice")
#'gene<-results(myDSResults)[,"gene"][1]
#'##Plot gene results filtering low expressed isoforms
#'g<-plotGeneResults(myDSResults, gene)
#'if(interactive()){
#'g
#'}
#'##Plot gene results keeping low expressed isoforms
#'g<-plotGeneResults(myDSResults, gene, filterLowExpIso=FALSE)
#'if(interactive()){
#'g
#'}
#'##Plot isoform bar plots keeping low expressed isoforms
#'g<-plotGeneResults(myDSResults, gene, filterLowExpIso=FALSE, group=FALSE)
#'if(interactive()){
#'g
#'}
setGeneric(name="plotGeneResults", def=function(myNBRes, gene,filterLowExpIso=
    TRUE, filterNotSignificant=TRUE,adjusted=TRUE, p.value=0.05, group=TRUE){
        standardGeneric("plotGeneResults")
})
#'@name plotGeneResults
#'@rdname NBSpliceRes-plotGeneResults
#'@aliases plotGeneResults,NBSpliceRes-method
#'@inheritParams plotGeneResults
setMethod(f="plotGeneResults", signature="NBSpliceRes", definition=function(
    myNBRes, gene, filterLowExpIso=TRUE, filterNotSignificant=TRUE,
    adjusted=TRUE, p.value=0.05,group=TRUE){
        DSDF<-GetGeneResults(myNBRes, gene, filterLowExpIso)
        DSDF$Filt<-!is.na(DSDF[, "pval"])
        if(filterNotSignificant){
            if(adjusted){
                DSDF$FiltDS<-DSDF[, "FDR"]< p.value
            }else{
                DSDF$FiltDS<-DSDF[, "pval"]< p.value
            }
            DSDF<-DSDF[DSDF$FiltDS,]
        
        }
        if(filterLowExpIso){
            
            DSDF<-DSDF[DSDF$Filt,]
        
        }
        
        rsh<-melt(DSDF[, seq_len(4)], id.vars=c("iso", "gene"))
        condVars<-colnames(DSDF)[grep("ratio", colnames(DSDF))]
        value<-variable<-iso<-NULL
        g<-ggplot(rsh, aes(x=variable, y=value,fill=iso))
        if(group){
            g<-g+geom_bar(stat="identity")          
        }else{
            g<-g+geom_bar(aes(group=iso),stat="identity", position="dodge")
        }
        g<-g+labs(x="Condition", y="Relative expression",fill="Isoform")
        g<-g+theme(panel.background=element_rect(fill="white", color="black"),
            legend.key=element_rect( fill="white", color="white"), 
            panel.grid.major=element_line(color="lightgoldenrod3"), 
            panel.grid.minor=element_line(color="tomato", linetype = "dashed"),
            legend.text = element_text(size = 12),legend.title = element_text(
            size = 14, face = "bold"), axis.text=element_text(size=12), 
            plot.title=element_text(size=14, hjust=0.5), axis.title=
            element_text(size=14))
        return(g)
})