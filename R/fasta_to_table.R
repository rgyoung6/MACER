#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'

# This function takes the two line or multiline fasta and puts it into a table
##################################### fasta_to_table FUNCTION ##############################################################
fasta_to_table <- function(Seq_file){


  #Initializing the flag
  fasta_flag=0
  Header<-NULL
  Sequence<-NULL

  for (j in 1:nrow(Seq_file)){

    #This is setting up a flag so the first time we initialize the matrix and then the second time we rbind to the matrix for the fasta file format
    if (fasta_flag==0){
      if(grepl(">",Seq_file[j,1])==TRUE){
        Header<-as.vector(Seq_file[j,1])
        seq_concate<-""
        fasta_flag=1
      }
    }else{
      if(grepl(">",Seq_file[j,1])==TRUE){
        Sequence<-rbind(Sequence,seq_concate)
        Header<-rbind(Header,as.vector(Seq_file[j,1]))
        seq_concate<-""
      }else{
        seq_concate<-paste(seq_concate,Seq_file[j,1],sep="")
      }
    }
  }
  Sequence<-rbind(Sequence,seq_concate)

  two_column_fasta<-cbind(Header,Sequence)

  return(two_column_fasta)

}

