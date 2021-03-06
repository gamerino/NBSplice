% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NBSpliceRes-initialize.R
\docType{methods}
\name{NBSpliceRes-initialize}
\alias{NBSpliceRes-initialize}
\alias{initialize,NBSpliceRes-method}
\title{NBSpliceRes object constructor.}
\usage{
\S4method{initialize}{NBSpliceRes}(.Object, results, lowExpIndex, contrast,
    dispersion)
}
\arguments{
\item{.Object}{NBSpliceRes class object.}

\item{results}{Data.frame with NBTest results of expressed isoforms.}

\item{lowExpIndex}{Numeric indicating the positions of low expressed isoforms.}

\item{contrast}{Character indicating the contrast used for NBTest.}

\item{dispersion}{Numeric with the estimated gene dispersions.}
}
\value{
NBSpliceRes object.
}
\description{
\code{initialize} creates an NBSpliceRes object
}
\note{
see full example in \code{\link{NBSpliceRes-class}}
}
\examples{

data(myDSResults, package="NBSplice")
myResults<-results(myDSResults, filter=FALSE)
myLowExpIdx<-lowExpIndex(myDSResults)
myContrast<-contrast(myDSResults)

myNewDSResults<-NBSpliceRes(myResults, myLowExpIdx, myContrast)

}
\seealso{
\code{\link{NBSpliceRes}}

Other NBSpliceRes: \code{\link{GetDSGenes}},
    \code{\link{GetDSResults}}, \code{\link{GetGeneResults}},
    \code{\link{NBSpliceRes-class}},
    \code{\link{NBSpliceRes}}, \code{\link{myDSResults}},
    \code{\link{plotGeneResults}},
    \code{\link{plotRatiosDisp}}, \code{\link{plotVolcano}}
}
\author{
Gabriela A. Merino \email{merino.gabriela33@gmail.com} and Elmer A. 
Fernandez \email{efernandez@bdmg.com.ar}
}
