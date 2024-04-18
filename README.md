# TMP

TMP leverages PCA and NNLS (Non-Negative Least Squares) to assess the ancestral composition of admixed individuals with high accuracy and reliability.

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
```

## Introduction

ASAP estimations are based on a PCA where the target and source groups are available. Along with those, we suggest to use additional groups to better define the PC space. Once that the PC space is defined, a set of NNLS is then applied on the PC coordinates, effectively summarizing the genetic ancestry. 

For its most basic usage ASAP needs:  

- A PCA matrix, a dataframe with N PCs
- A list of target and source groups (or samples)

### Basic Usage Example
  
```{r, eval=FALSE}
$ pca = read_eigen(pca_input = 'data/TOY.pca.evec')
$ ASAP_result = ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
```

ASAP package has two functions to read the PCA matrix, ```read_eigen()``` and ```read_flash()```.  
* ```read_eigen()``` will read a PCA that has been created with smartpca from the [EIGENSOFT](https://github.com/DReichLab/EIG) software.  
* ```read_flash()``` instead, will read a PCA that has been created with [flashpca](https://github.com/gabraham/flashpca) software.  
The goal of both functions is to set the PCA file as follows:  
  
| POP | IND   | PCN |
| --- | ----  | ------- |
| ASW | ASW_1 | 0.01 |
| ASW | ASW_2 | 0.02 |
| BEB | BEB_1 | 0.08 |
  
If neither ```read_eigen()``` nor ```read_flash()``` is for you, you might want to simply use ```read.table()```, and set the file so that it has the aforementioned look.  

For the sake of the example, let's say you obtained a PCA from the software EIGENSOFT, using smartpca.  

```{r, eval=FALSE}
$ pca = read_eigen(pca_input = 'data/TOY.pca.evec')
```  
  
The funtion ```ASAP()``` requires also a list of the target and reference groups: you can provide the list in two ways.

1) You can provide an 'AS_file': a file with the list of the Admixed groups (A) and the Source groups (S). If you want to use ASAP sample-wise rather than group-wise, simply adjust the PCA file so that in the 'POP' column is identical to the 'IND' column, and set the AS_file with the samples list, rather than the group list.   
The AS_file is a two-columns file with the population list on the first column, and the 'A/S' information on the second column. The 'A/S' information stands for Admixed (A) or Source (S). For each population/group we will indicate whether ASAP should consider it as a Source (S) or as an admixed target (A), the file looks like this:  

| POP | A/S |
| --- | --- |
| ESN | S |
| CEU | S |
| GIH | A |

To read the AS_file, a simple ```read_table(file, header=T)``` will be sufficient.   

```{r, eval=FALSE}
$ AS_file = read.table('data/Example_AS_eigen', header=TRUE)
```

With the PCA and AS_file loaded, we are finally ready to run ASAP as follows:

```{r, eval=FALSE}
$ ASAP_result = ASAP(pca_input = pca, as_file = AS_file)
```
  
2) You can avoid relying on the AS_file if you wish, using a vector of the target and source groups directly in ASAP() function, as follows:
  
```{r, eval=FALSE}
$ ASAP_result = ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
```
  
Finally, if you want to save ASAP results on a table-like format, you can use ```write_ASAP()```, this way:
  
```{r, eval=FALSE}
$ ASAP_result <- ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
$ write_ASAP(ASAP_input = ASAP_result, output_name = 'my_dir/my_asap_results.txt')
```

## Cite TMP

If you use TMP, please cite [BiorXiv preprint]()

## Contact

For questions and bug reports please contact [LM](mailto:ludovica.molinaro@kuleuven.be).
