Installation requires some additional packages.  The following code should install everything required:

```{r}
install.packages(c('BiocManager', 'remotes'))
BiocManager::install(c("BiocStyle","GGally","dplyr","ggplot2","data.table","lmerTest",
           "Biostrings","Rsubread","ShortRead","gridExtra","tibble","tidyr",
           "optimx","stringr","hadley/multidplyr"))
BiocManager::install("grimbough/MEMcrispR")           
```