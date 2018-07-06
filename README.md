# NBSplice
NBSplice
=======

The `NBSplice` R package allows the evaluation of differential gene splicing through negative binomial generalized linear models fitted at the gene level using isoform expression counts.
To install and load the `NBSplice` package in R:

```r
if (!require("BiocManager"))
install.packages("BiocManager")
BiocManager::install("NBSplice", version="devel")
library("NBSplice")
```

To find the vignette associated with `NBSplice` using R:

```r
browseVignettes(package = "NBSplice")
```

