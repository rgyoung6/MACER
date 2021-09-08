#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'

# Identify the forward and reverse beginning of the reference sequence in relation to the output alignment
############################## forward_reverse_trim FUNCTION ############################
forward_reverse_trim <- function(seq_data_matrix){

  forward_seq_start = NULL
  reverse_seq_start = NULL

  #Identifying the beginning of the sequence alignment by identifying the column with non-gaps using per_align_start
  for (i in 3:ncol(seq_data_matrix)){
    if(is.null(forward_seq_start)){
      if(nrow(subset(seq_data_matrix, (seq_data_matrix[,i]=="-")|(seq_data_matrix[,i]=="N")))==0 ){
        forward_seq_start = i
        break;
      }
    }
  }

  #Identifying the end of the sequence alignment by identifying the column with 100% non-gaps
  for (i in ncol(seq_data_matrix):3){
    if(is.null(reverse_seq_start)){
      if(nrow(subset(seq_data_matrix, (seq_data_matrix[,i]=="-")|(seq_data_matrix[,i]=="N")))==0 ){
        reverse_seq_start = i
        break;
      }
    }
  }

  forward_reverse_pos<-c(forward_seq_start,reverse_seq_start)

  return(forward_reverse_pos)

}
