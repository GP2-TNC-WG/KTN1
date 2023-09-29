## Power calculations for expression data

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RNASeqPower")
library("RNASeqPower")

#0 Months
rnapower(depth=25, n=1568, n2=1042, cv=0.4, effect=2, alpha=.05)

#24 Months
rnapower(depth=25, n=785, n2=328, cv=0.4, effect=2, alpha=.05)