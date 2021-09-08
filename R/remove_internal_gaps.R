#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'

# Identify and remove internal gaps
############################## remove_internal_gaps FUNCTION ############################
remove_internal_gaps <- function(output_fasta_matrix, pigl){

  #initializing the variable
  col_to_remove = NULL

  #finding columns with the number sequences with a gap at that position greater than the internal gap percentage, then remove the sequences with
  # nucleotides causing that internal gap and then add the column to a vector for removal as it would now only have gaps
  for(k in 3:ncol(output_fasta_matrix)){
    if(nrow(subset(output_fasta_matrix, (output_fasta_matrix[,k]=="-")|(output_fasta_matrix[,k]=="N")))>(nrow(output_fasta_matrix)*pigl)){
      #remove the sequence(s) with the nucleotide
      output_fasta_matrix<-subset(output_fasta_matrix, (output_fasta_matrix[,k]=="-")|(output_fasta_matrix[,k]=="N"))
      #setting up a vector to hold all of the rows to be removed note that I am not subtracting 1 as I need to account for the row with the header
      col_to_remove<-c(col_to_remove,(k))
    }
  }
  if (!is.null(col_to_remove)){
    #deleting the internal gap sequences
    output_fasta_matrix<-output_fasta_matrix[,-c(col_to_remove)]
  }

  return(output_fasta_matrix)

}
