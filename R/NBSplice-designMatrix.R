#'An example of a design matrix helpful to demonstrate
#'the use of the NBSplice R package.
#'
#'A data.frame representing the design matrix related to the isoCounts dataset.
#'
#'\describe{
#'This data.frame complements the expression matrix provided as the
#'\link{isoCounts} data. See \link{isoCounts} man page for information about
#'\link{designMatrix} construction.
#'Each row of this data.frame specifies one experiment sample.
#'\item{samples}{Column with samples names}.
#'\item{condition}{Column indicating the condition of each sample}.
#'}
#'
#'@include NBSplice-geneIso.R
#'@docType data
#'@format A data.frame object
#'@source see \code{\link{IsoDataSet-class}}
#'@name designMatrix
#'@family IsoDataSet
#'@author Gabriela A. Merino \email{merino.gabriela33@@gmail.com} and 
#'Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
NULL