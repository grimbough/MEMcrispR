Installation requires some additional packages.  The following code should install everything required:

```{r}
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("BiocStyle","GGally","dplyr","ggplot2","data.table","lmerTest",
           "Biostrings","Rsubread","ShortRead","gridExtra","tibble","tidyr",
           "tools","optimx","stringr","hadley/multidplyr"))
biocLite("grimbough/METAcrispR")           
```