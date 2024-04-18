#' asap
#'
#' asap function takes a PCA as input and analyses it via NNLS to describe admixed individuals as a mixture of sources groups.
#' @param pca_input R data.frame of PCA with N PCs
#' @param as_file R data.frame with two columns: POP and A/S, POP column lists all populations to be considered, A/S indicatea whether the population should be considered as Admixed ('A') or as Source ('S')
#' @param sources R vector indicating the groups that should be considered as Sources
#' @param admixed R vector indicating the groups that should be considered as Admixed
#' @examples
#' pca = read.table('example.pca') #OR
#' pca = read_eigen('data/TOY.pca.evec') #OR
#' pca = read_flash('data/TOY_flash.pca')
#'
#' example_AS = read.table('data/Example_AS', header = T)
#'
#' asap(pca_input = pca, as_file = example_AS) #OR
#' asap(pca_input = pca, sources = c('Source1','Source2','Source3'), admixed = c('Admixed1','Admixed2')
#' @return Returns a list containing the ancestries proportions per each Admixed group
#' @export

asap <- function(pca_input, as_file, sources = NULL, admixed = NULL) {

  if(missing(pca_input)) {
    stop('ERROR: asap() needs a PCA matrix as input')
  }

  else if(is.null(sources)&is.null(admixed)&missing(as_file)) {

    stop('ERROR: asap() needs a reference file to assign sources and admixed samples. Assign a file to as_file, OR a vector to source_list AND vector to admixed_list arguments')

  } else if(is.null(sources)&is.null(admixed)&!missing(as_file)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and average their PCs values

    ## Aggregate
    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'

    ## Group_by
    #pca_grouped_avg <- subset_pca %>% group_by(POP) %>% summarise_all(funs(mean=mean))
    #pca_grouped_int = as.matrix(sapply(pca_grouped_avg, as.numeric))

    ## Selecting Sources and Admixed based on as_file file

    Admixed_pop = as.character(as_file$POP[as_file$A.S == 'A'])
    Source_pop = as.character(as_file$POP[as_file$A.S == 'S'])

    ### Source copying vector

    Source_cv = pca_aggregated[pca_aggregated$POP %in% Source_pop,]
    if (nrow(Source_cv) == 0) {stop("Your source vector is empty, check your input files")}
    if(any(is.na(Source_cv))) {
      print(Source_cv)
      stop("The source vector contains NA, check your input files")}

    ### Admixed copying vector

    Admixed_cv = pca_aggregated[pca_aggregated$POP %in% Admixed_pop,]
    if (nrow(Admixed_cv) == 0) {stop("Your admixed vector is empty, check your input files")}
    if(any(is.na(Admixed_cv))) {
      print(Admixed_cv)
      stop("The admixed vector contains NA, check your input files")}

    ## NNLS analyses
    asap_RES = nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cold and rows for nnls.mat
    #colnames(asap_RES) = Source_pop
    #rownames(asap_RES) = Admixed_pop

    ## Rename cols and rows for nnls.mat2
    colnames(asap_RES[[1]]) = Source_cv$POP
    rownames(asap_RES[[1]]) = Admixed_cv$POP

    ## printing/returning output matrix
    return(asap_RES)
  }  else if(missing(as_file)&!is.null(sources)&!is.null(admixed)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and average their PCs values

    ## Aggregate
    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'

    ## Group_by
    #pca_grouped_avg <- subset_pca %>% group_by(POP) %>% summarise_all(funs(mean=mean))
    #pca_grouped_int = as.matrix(sapply(pca_grouped_avg, as.numeric))

    ## Selecting Sources and Admixed based on input lists

 ### Source copying vector

    Source_cv = pca_aggregated[pca_aggregated$POP %in% sources,]
    if (nrow(Source_cv) == 0) {stop("Your source vector is empty, check your input files")}
    if(any(is.na(Source_cv))) {
      print(Source_cv)
      stop("The source vector contains NA, check your input files")}

  ### Admixed copying vector
    Admixed_cv = pca_aggregated[pca_aggregated$POP %in% admixed,]
    if (nrow(Admixed_cv) == 0) {stop("Your admixed vector is empty, check your input files")}
    if(any(is.na(Admixed_cv))) {
      print(Admixed_cv)
      stop("The admixed vector contains NA, check your input files")}

    ## NNLS analyses
    asap_RES = nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cols and rows
    colnames(asap_RES[[1]]) = Source_cv$POP
    rownames(asap_RES[[1]]) = Admixed_cv$POP

    ## printing/returning output matrix
    return(asap_RES)
  } else {
    stop('ERROR: optional arguments were not found.')
  }
}


#' write_asap
#'
#' write_asap allows to save ASAP results in a table-like format.
#' @param asap_input R list returned by asap() function
#' @param output_name string containing the file output name
#' @examples
#' asap_results <- asap(pca_input = pca, as_file = example_as)
#' write_asap(asap_input = asap_results, output_name = 'my_dir/my_asap_results.txt')
#' @export

write_asap <- function(asap_input,output_name) {

  df_asap = data.frame(asap_input[[1]])
  df_residuals = data.frame(asap_input[[2]])
  names(df_residuals) = 'Residuals'

  final_df = cbind(df_asap,df_residuals)

  write.table(final_df,file=output_name)
}
