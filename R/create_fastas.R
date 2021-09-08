#Roxygen2 Documentation:

#' @export
#'
#' @title Table To FASTA
#'
#' @author Rekkab Singh Gill
#'
#' @description Using the output table from the download script and the user built genus-marker name parameter file to take the downloaded data and place them into fasta files.
#'
#' @details
#' Input: File with list of genera with the molecular markers names below the taxa. The information to create this parameters file can be obtained from A_Summary.txt file from the download script results.
#' For further details please see the documentation.
#'
#' @examples
#' \dontrun{
#' create_fastas()
#' create_fastas(no_marker = TRUE, no_taxa = TRUE)
#' create_fastas(no_seq  = TRUE, name_issue = TRUE)
#' }
#'
#' @param data_file NULL prompts the user to indicate the location of the data file in the format of the auto_seq_download output, anything other than NULL then the string supplied will be used for the location; default NULL
#' @param input_file NULL prompts the user to indicate the location of the input file used to select through point and click prompts, anything other than NULL then the string supplied will be used for the location; default NULL
#' @param output_folder NULL prompts the user to indicate the location of the output file through point and click prompts, anything other than NULL then the string supplied will be used for the location; default NULL
#'
#' @param no_marker If set to TRUE then will include records filtered out due to no marker data. Default is FALSE to not include records with no marker data.
#' @param no_taxa If set to TRUE then will include records filtered out due to no taxa data. Default is FALSE to not include records with no taxa data.
#' @param no_seq If set to TRUE then will include records filtered out due to no sequence data. Default is FALSE to not include records with no sequence data.
#' @param name_issue If set to TRUE then will include records filtered out due to genus and species names with more than two terms. Default is FALSE to not include records with taxonomic naming issues.
#' @param taxa_digits If set to TRUE then will include records filtered out due to genus or species names containing digits. Default is FALSE to not include records with digits in the taxonomic naming.
#' @param taxa_punct If set to TRUE then will include records filtered out due to the presence of punctuation in the genus or species names. Default is FALSE to not include records with punctuation in the taxonomic naming.
#' @returns
#' This script outputs a fasta file of sequences for each column in the submitted parameters file. These files are named with the genera of interest and the first marker name in the column of the parameters file.
#' These files are located in the folder where the Total_tables.txt file is located.
#'
#' @references
#' <https://github.com/rgyoung6/MACER>
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'
#' @seealso
#' create_fastas()
#' align_to_ref()
#' barcode_clean()
#'
#' @import utils
#'


#********************************************Main program section***********************************************
##################################### Main FUNCTION ##############################################################
create_fastas <- function(data_file = NULL, input_file = NULL, output_folder = NULL, no_marker = FALSE, no_taxa = FALSE, no_seq= FALSE, name_issue = FALSE, taxa_digits = FALSE, taxa_punct = FALSE ){


  #Check to see if a path to the data_file of interest was submitted in the function call
  if (is.null(data_file)){

    n <- substr(readline(prompt="Please select the total tables file.  Hit enter key to continue..."),1,1)
    #prompting the user for the file through file.choose
    total_tables_file<-file.choose()

  }else {

    total_tables_file=data_file

  }

  if (is.null(input_file)){

    n <- substr(readline(prompt="Please select the file with genus and the list of molecular markers of interest. Hit enter key to continue..."),1,1)
    #prompting the user for the file through file.choose
    genus_marker_file<-file.choose()

  }else{

    genus_marker_file=input_file

  }

  if (is.null(output_folder)){

    current_path <- dirname(total_tables_file)

  }else{

    current_path=output_folder

  }

  #read in the data from the total tables file
  total_table_data<-as.data.frame(read.table(total_tables_file,header=TRUE,sep="\t",na.strings="" ))

  #check to see its not empty
  if(nrow(total_table_data) < 1)
  {
    stop('Unable to create fasta file, input file not valid')
  }

  #Get the columns of interest to collapse into a header
  carrot<-as.data.frame(rep(">", nrow(total_table_data)))

  #Collapse the header table into a header vector
  headers<-as.data.frame(apply(total_table_data[,c(3,4,6,7,8,9)],1,paste,collapse="|"))

  #Create a fasta header and place it at the end of table
  headers<-cbind(carrot, headers)

  #combine the headers into a single column
  headers<-as.data.frame(paste0(headers[,1],headers[,2]))

  #Change the column name
  colnames(headers) <- c("Header")

  #Add the headers on the end of the total_table_data
  total_table_data<-cbind(total_table_data,headers)

  #read in the data from the total tables file
  genus_marker_data<-as.data.frame(read.table(genus_marker_file,header=FALSE,sep="\t",na.strings="" ))

  #check to see its not empty
  if(nrow(genus_marker_data) < 1)
  {
    stop('Unable to create fasta file, genus and marker names file not valid')
  }

  #see what parameters are set by the user and add them to the parameter list
  params_list <- list('-')
  if(no_marker == TRUE) params_list <- append(params_list,'No_Marker')
  if(no_taxa == TRUE) params_list <- append(params_list,'No_Taxa')
  if(no_seq == TRUE) params_list <- append(params_list,'No_Seq')
  if(name_issue == TRUE) params_list <- append(params_list,'Taxa_name_issue')
  if(taxa_digits == TRUE)params_list <- append(params_list,'Taxa_w_digits')
  if(taxa_punct == TRUE) params_list <- append(params_list,'Taxa_w_punct')

  #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
  total_table_data_out<-total_table_data[(total_table_data$Flags %in% params_list),]

  #check to see its not empty
  if(nrow(total_table_data_out) < 1)
  {
    stop('Unable to create fasta file, no records that meet the parameters submitted')
  }

  #Getting all of the column names to loop through these genera
  genus_col_names<-genus_marker_data[1,]

  #Removing the first row off the genus_col_names table
  genus_marker_data<-as.data.frame(genus_marker_data[-1,])

  #Using the genus and marker file create fastas for each column of interest
  for(genus_col in 1:length(genus_col_names)){

    # genus_col=1

    #Get the records with the genus of interest
    genus_data_out<-total_table_data_out[(total_table_data_out$Genus %in% genus_col_names[genus_col]),]

    #check to see its not empty
    if(nrow(total_table_data_out) < 1)
    {
      print(paste0(genus_col_names[genus_col], " - Unable to create fasta file, no records for this genus"))

    }else{

      #get the molecular marker names of interest
      marker_names<-unlist(genus_marker_data[,genus_col])

      #Remove the NA for the molecular marker names
      marker_names<-marker_names[!is.na(marker_names)]

      #Get the records for the molecular markers of interest
      genus_data_out<-genus_data_out[(genus_data_out$Gene %in% marker_names),]

      #check to see its not empty
      if(nrow(genus_data_out) < 1)
      {
        print(paste0(genus_col_names[genus_col], " - Unable to create fasta file, no records for this genus and molecular markers"))
      }else{

        #select only the header and the sequence columns
        genus_data_out<-cbind(genus_data_out$Header,genus_data_out$Sequence)

        #write out the final_fasta_frame to a fasta file
        newFilePath <- file.path(current_path, paste0(genus_col_names[genus_col], "_", marker_names[1]  ,".fas"))
        write.table(genus_data_out,newFilePath, na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=FALSE)

        (print("Complete, please see the file location with your input table for results."))

      }#End of if statement checking if there are molecular marker data for the genus of interest

    }#End of if statement checking to see if the data are empty at the genus level

  }#End fo the for loop

}# End of the function
