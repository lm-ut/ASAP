## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ASAP)

## ---- eval=FALSE--------------------------------------------------------------
#  $ pca = read_eigen(pca_input = 'data/TOY.pca.evec')
#  $ ASAP_result = ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))

## ---- eval=FALSE--------------------------------------------------------------
#  $ pca = read_eigen(pca_input = 'data/TOY.pca.evec')

## ---- eval=FALSE--------------------------------------------------------------
#  $ AS_file = read.table('data/Example_AS_eigen', header=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  $ ASAP_result = ASAP(pca_input = pca, as_file = AS_file)

## ---- eval=FALSE--------------------------------------------------------------
#  $ ASAP_result = ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))

## ---- eval=FALSE--------------------------------------------------------------
#  $ ASAP_result <- ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
#  $ write_ASAP(ASAP_input = ASAP_result, output_name = 'my_dir/my_asap_results.txt')

## ---- eval=FALSE--------------------------------------------------------------
#  $ ASAP_result <- ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
#  $ plot_asap(ASAP_result, 'ASAP_plot', type_ggplot = 'YES')

## ---- eval=FALSE--------------------------------------------------------------
#  $ ASAP_result <- ASAP(pca_input = pca, sources = c('ESN','CEU'), admixed=c('GIH'))
#  $ plot_asap(ASAP_result, 'ASAP_plot')

## ---- eval=FALSE--------------------------------------------------------------
#  $ source_pairs = read.table('see_example_in_II', header = T)
#  $ pcs_distances(pca_input = pca, output_name = 'ED_cumuldistances', sources_file = source_pairs, return_plot = 'YES')

## ---- eval=FALSE--------------------------------------------------------------
#  $ pc_target <- c()
#  $ for(i in 1:23) {pc_target[[i]] <- paste0("PC",i)}
#  $ cols_target = c("POP","IND",pc_target)
#  $ ASAP_result <- ASAP(pca_input = pca[,cols_target], sources = c('Source1','Source2'), admixed=c('Admixed1'))

## ---- eval=FALSE--------------------------------------------------------------
#  $ AS_file = read.table('data/Example_AS_eigen', header=TRUE)
#  $ ASAP_resampling = read.resampling(path_tofiles = 'where/I/stored/my/PCAs/', file_pattern = 'PCA_Resampled_', as_file = AS_file, eigentype)

## ---- eval=FALSE--------------------------------------------------------------
#  $ main_pca = read_eigen(pca_input = 'data/TOY.pca.evec')
#  $ AS_file = read.table('data/Example_AS_eigen', header=TRUE)
#  $ ASAP_main_result <- ASAP(pca_input = main_pca, as_file = AS_file)
#  $ ASAP_resampling = read.resampling(path_tofiles = 'where/I/stored/my/resampled/PCAs/', file_pattern = 'PCA_Resampling_', as_file = AS_file, eigentype)
#  $ chromovec = 20
#  $ ASAP_SE = se.resampling(ASAP_main_result, ASAP_resampling, chromovec)

