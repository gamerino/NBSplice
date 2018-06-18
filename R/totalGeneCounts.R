#'@title
#'Auxiliary function to build gene expression matrix
#'@description
#'\code{totalGeneCounts} obtains the total counts for a gene in a sample by the
#'sum of the counts of all its isoforms. This count is performed before any
#'filter step to avoid erroneus estimation of isoform ratios.
#'
#'@param isoCounts Matrix containing the isoform expression values with
#'isoforms in rows and samples in columns. Object row names should 
#'be the names of the isoforms and the column names, the names of each sample.
#'@param geneIso Data.frame with two columns called "isoform_id" and "gene_id"
#'especifying which are the isoform (rows of the expression matrix) for each
#'gene
#'@param BPPARAM A bpparam object specifying the number of cpus to be used.
#'
#'@return A matrix of the same size as the isoCounts matrix with the total 
#'counts for each gene in each sample in CPM scale.
#'
#'@include IsoDataSet.R
#'@export totalGeneCounts
#'@docType methods
#'@name totalGeneCounts
#'@rdname totalGeneCounts
#'@aliases totalGeneCounts-methods
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and Elmer A.
#'Fernandez \email{efernandez@bdmg.com.ar}
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bplapply
#'@importFrom BiocParallel bpparam
#'@import edgeR
#'@importFrom edgeR cpm
#'
#'@examples
#'
#'data(isoCounts, package="NBSplice")
#'data(geneIso, package="NBSplice")
#'
#'myGeneCounts<-totalGeneCounts(isoCounts, geneIso)
#'
totalGeneCounts<-function(isoCounts, geneIso, BPPARAM=bpparam()){
    # CPM normalization
    isoCounts<-cpm(isoCounts)
    genes<-as.character(unique(geneIso[,"gene_id"]))

    allCounts<-do.call(rbind,bplapply(seq_along(genes), 
    function(x){
        gene_cm<-isoCounts[as.character(geneIso[,"isoform_id"][geneIso[,
            "gene_id"] == genes[x]]), , drop=FALSE]  
        total<-matrix(1, nrow=nrow(gene_cm), ncol=1) %*% matrix(
            colSums(gene_cm), nrow=1)
        dimnames(total)<-dimnames(gene_cm)
        return(total)
    }, BPPARAM=BPPARAM))
    colnames(allCounts)<-paste(colnames(isoCounts), "_All", sep="")
    return(allCounts)
}