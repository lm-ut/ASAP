# ASAP

ASAP leverages PCA and NNLS (Non-Negative Least Squares) to assess the ancestral composition of admixed individuals with high accuracy and reliability.

## Installation

```{r, eval=FALSE}
install.packages("devtools") 
devtools::install_github("lm-ut/tmp", dependencies = TRUE)
library("tmp")
```

ASAP requires the following R packages: 

```{r, eval=FALSE}
install.packages("dplyr")
install.packages("nnls")
install.packages("ggplot2")
