% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NBSpliceRes-print.R
\docType{methods}
\name{NBSpliceRes-print}
\alias{NBSpliceRes-print}
\alias{print,NBSpliceRes-method}
\title{Print an NBSpliceRes object.}
\usage{
\S4method{print}{NBSpliceRes}(x, ...)
}
\arguments{
\item{x}{NBSpliceRes class object.}

\item{...}{Included for generic print compatibility.}
}
\value{
Console output of the object.
}
\description{
Generic \code{print} method for NBSpliceRes class and descendants.
}
\note{
see full example in \code{\link{NBSpliceRes-class}}
}
\examples{

data(myDSResults, package="NBSplice")

print(myDSResults)

}
\author{
Gabriela A. Merino \email{merino.gabriela33@gmail.com},
and Elmer A. Fernandez \email{efernandez@bdmg.com.ar}
}
