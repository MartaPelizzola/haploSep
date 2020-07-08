# haploSep

The R package haploSep provides functions to estimate multiple haplotypes from allele frequency data. In order to use this package allele frequencies from several samples are needed. With the function haploSelect the number of dominant haplotypes in a population can be estimated. The function haploSep provides estimates for the haplotype structure and for their frequencies for each sample, together with the corresponding quality scores. 

haploSep also provides a customized plot function and a function to compare different haplotype data sets.

For more information about haploSep see our manuscript on bioRxiv "Multiple Haplotype Reconstruction from Allele Frequency Data" by Pelizzola et al.  

# Installation 

Copying the following line in R one should be able to install and load the package:

```{r}
install.packages('devtools')
library(devtools)
devtools::install_github("MartaPelizzola/haploSep")

library(haploSep)
```

# Authors
Merle Behr (merle.behr@berkeley.edu)
Housen Li (housen.li@mathematik.uni-goettingen.de)
Marta Pelizzola (marta.pelizzola@vetmeduni.ac.at)
