#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'
#' @import utils
#'
############################################ BOLD DOWNLOAD FUNCTION ####################################################
BOLD_download <- function(target_genera, BOLD_folder_str, main_file_folder)
{

  #bring in the genera of interest and the location of the filefolder where the output will be stored
  #******************************************** BOLD API DOWNLOAD ***********************************************
  #For species names get the BOLD full record file for all molecular markers
  html_string<-paste0("http://v4.boldsystems.org/index.php/API_Public/combined?taxon=", target_genera, "&format=tsv")
  BOLD_out_name<-paste0(BOLD_folder_str,"/",target_genera,"_BOLD_download.txt")
  download.file(html_string,BOLD_out_name, quiet=TRUE)
  #************************************************************ PRODUCE TABLE FROM BOLD DATA DOWNLOAD ***********************************

  #Getting the information about the downloaded file to check if there is any data downloaded.
  BOLD_download_file_info<-file.info(BOLD_out_name)

  #Here the if statement checks to see if there is data in the file info vector. The first column represents the size fo the file of interest
  if(BOLD_download_file_info[1]!=0){

    # load in the data file
    BOLD_download_data_table<-read.table(BOLD_out_name,header=T,sep="\t",dec=".", quote="",comment.char = "", fill =TRUE)

    # Create a vector with BOLD for the number of rows of the table
    DB<-c(rep("BOLD",nrow(BOLD_download_data_table)))

    # creating the output standard data table to combine with GenBank
    BOLD_download_data_table<-cbind(DB,BOLD_download_data_table[,c(1,71,22,8,70,72)])

    #remove the BOLD: from the BIN column
    BOLD_download_data_table[,5]<-gsub("BOLD:","",BOLD_download_data_table[,5])

  } else {

    BOLD_download_data_table=NULL

  } #Closing the if-else statement that checked to see if there were records downloaded

  return(BOLD_download_data_table)

}
