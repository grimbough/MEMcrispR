Installation requires some additional packages.  The following code should install everything required:

```{r}
install.packages(c('BiocManager', 'remotes'))
BiocManager::install("grimbough/MEMcrispR", dependencies = TRUE)           
```
