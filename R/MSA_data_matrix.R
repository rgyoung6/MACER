#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'

# This function takes the data table (header, sequence) and make it into a data matrix with each nucleotide position of the MSA in a column
##################################### MSA_data_matrix FUNCTION ##############################################################
MSA_data_matrix <- function(seq_data_table){

  #grabbing the headers and placing them in the outgoing matrix
  seq_data_matrix <- seq_data_table[,1, drop=FALSE]

  # The t is transposing the matrix. strsplit the second column containing the sequence data in to columns for each character and bind to the headers in the outgoing matrix seq_data_file
  v <- t(as.data.frame(strsplit(as.character(seq_data_table[,2]), "")))
  rownames(v)<-NULL
  seq_data_matrix<-cbind(seq_data_matrix,as.character(seq_data_table[,2]),v)

  #returning the final table as output from the function
  return(seq_data_matrix)

}
