#Roxygen2 Documentation:

#' @export
#'
#' @title Automatic Sequence Download
#'
#' @author Robert G. Young
#'
#' @description
#' Takes a list of genera, as supplied by the user, and searches and downloads molecular sequence data from BOLD and Genbank.
#'
#' @details
#' User Input: A list of genera in a text file in a single column with a new line at the end of the list.
#'
#' @examples
#' \dontrun{
#' auto_seq_download()
#' auto_seq_download(BOLD_database = TRUE, NCBI_database = FALSE)
#' auto_seq_download(BOLD_database = FALSE, NCBI_database = TRUE)
#' }
#'
#' @param BOLD_database TRUE is to include, FALSE is to exclude; default TRUE
#' @param NCBI_database TRUE is to include, FALSE is to exclude; default TRUE
#' @param search_str NULL uses the default string, anything other than NULL then that string will be used for the GenBank search; default NULL.
#' The Default String is: (genus[ORGN]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])
#' @param input_file NULL prompts the user to indicate the location of the input file through point and click prompts, anything other than NULL then the string supplied will be used for the location; default NULL
#' @param output_file NULL prompts the user to indicate the location of the output file through point and click prompts, anything other than NULL then the string supplied will be used for the location; default NULL
#' @param seq_min holds the minimum length value to not flag the sequence; default 100
#' @param seq_max holds the maximum length value to not flag the sequence; default 2500
#'
#'
#' @returns Outputs: One main folder containing three other folders.
#' Main folder - Seq_auto_dl_TTTTTT_MMM_DD
#' Three subfolders:
#' 1. BOLD - Contains a file for every genus downloaded with the raw data from the BOLD system.
#' 2. NCBI - Contains a file for every genus downloaded with the raw data from GenBank.
#' 3. Total_tables - Contains files for the running of the function which include...
#' A_Summary.txt - This file contains information about the downloads.
#' A_Total_Table.tsv - A file with a single table containing the accumulated data for all genera searched.
#'
#' @references
#' <https://github.com/rgyoung6/MACER>
#' Young, R. G., Gill, R., Gillis, D., Hanner, R. H. (Submitted June 2021). Molecular Acquisition, Cleaning, and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal.
#'
#' @note
#' When using a custom search string for NCBI only a single genus at a time can be used.
#'
#' @seealso
#' create_fastas()
#' align_to_ref()
#' barcode_clean()
#'
#' @import httr
#' @import utils
#'

############################################ MAIN FOR FOR BOLD/NCBI DOWNLOAD ###################################################################################
auto_seq_download <- function(BOLD_database=TRUE, NCBI_database=TRUE, search_str= NULL, input_file= NULL, output_file=NULL, seq_min=100, seq_max=2500)
{

  #The following if checks to see if the user inputted a file location when calling the program if not then prompts the user to select
  if (is.null(input_file)){

    # prompting to choose the file of interest with the tab delimited info
    n <- substr(readline(prompt="Choose the file with the genera of interest to download. Hit enter key to continue..."),1,1)
    genera_list<-file.choose()

  } else if (!is.null(input_file)){

    genera_list=input_file

  }

  #assign the output if the user didn't supply an output location
  if (is.null(output_file)){

    Work_loc <- dirname(genera_list)
    print(genera_list)

  }else{

    Work_loc=output_file

  }

  # load in the data file in to the data for the lines
  genera_list<-read.table(genera_list,header=F,sep="\t",dec=".", fill = TRUE )

  #Check to see if there is a one column table with only genera names or a two column table including specific search criteria
  if(ncol(genera_list)==1){

    genera_list<-genera_list[,1]

  }else if (ncol(genera_list)==2){

    search_str<-genera_list[,2]
    genera_list<-genera_list[,1]

  }

  #Make a file folder to contain all results from the program
  main_file_folder<-paste0(Work_loc, "/Seq_auto_dl",format(Sys.time(), "_%H%M%S_%b_%d"))
  dir.create(main_file_folder)

  #Flags for the creation file folders to hold BOLD and GenBank outputs
  BOLD_flag = 1
  NCBI_flag = 1

  #Create a subfolder to keep all the genera tables
  table_folder<-paste0(main_file_folder, "/Total_Tables")
  dir.create(table_folder)

  #Initialize final output variables
  BOLD_data_table_cleaned = NULL
  NCBI_data_table_cleaned = NULL

  #Initialize the total_data_table variable that will be used for the final data frame
  total_data_table<- data.frame(uniqueID=character(), DB=character(), ID=character(), Accession=character(),Genus_Species=character(), Genus=character(),Species=character(),BIN_OTU=character(),Gene=character(),Sequence=character(),Flags=character(), stringsAsFactors=FALSE)

  #Initialize the output file for all of the molecular markers
  write.table(total_data_table ,file=paste0(table_folder,"/A_Total_Table.tsv"), na="", row.names=FALSE, col.names=TRUE, quote = FALSE,sep="\t", append=FALSE)

  #Initializing the file for output to the summary log file
  summary_log<-NULL

  #Initialize the output file for the summary log
  write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\t", append=TRUE)

  #Starting the initial loop to go through each of the items to search one at a time
  for (genera_list_loop_counter in 1:length(genera_list)) {

    summary_log<-paste0("Starting time - ", Sys.time(), " - ",genera_list[genera_list_loop_counter], " search ", genera_list_loop_counter, " of ",length(genera_list) )
    (print(summary_log))
    #Adding to the summary log for the running of the program and the results
    write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

    #This is a data table to hold the accumulation of all records from the target taxa list to use later to create marker specific fasta for all taxa
    total_data_table<-data.frame(DB=character(), ID=character(), Accession=character(), Genus=character(), Species=character(), BIN_OTU=character(), Gene=character(), Sequence=character(), stringsAsFactors=FALSE)

    #*****************************downloading from BOLD, attempt at least 3 times*******************************************************************************#

    if (BOLD_database==TRUE){

      #Adding the attempt time to the log file
      summary_log<-paste0("Attempting BOLD Download - ", Sys.time(), " - ", genera_list[genera_list_loop_counter])
      (print(paste0("Attempting BOLD Download - ", Sys.time(), " - ", genera_list[genera_list_loop_counter])))
      #Adding to the summary log for the running of the program and the results
      write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

      if (  BOLD_flag == 1){

        #Create a subfolder to keep all of the full record downloads from BOLD
        BOLD_folder_str<-paste0(main_file_folder, "/BOLD")
        dir.create(BOLD_folder_str)

        BOLD_flag = 0
      }
      #initialize the attempt variables, set the error occurrence flag to false and initialize an empty data frame
      attempt<- 1
      BOLD_data_table <- data.frame(matrix(nrow=0, ncol=2))

      while(attempt < 4)
      {
        if(attempt == 2)
        {
          # This is setting the internet wait time to a higher value than the 60 sec default. We needed to make it at least 120 as we were having timeout errors.
          options(timeout = 120)
        }
        else if (attempt == 3)
        {
          options(timeout = 240)
        }

        BOLD_data_table<-tryCatch({
          #Call the BOLD_download function send to it the location for the output file
          result<-as.data.frame(BOLD_download(genera_list[genera_list_loop_counter], BOLD_folder_str, main_file_folder))
        },
        error = function(e){
          print('Timeout Error, retrying with greater timeout limit')
        })

        #exit loop if error didn't occur
        if(typeof(BOLD_data_table) != 'character')
        {
          break
        }
        #increment the counter and reset the dataframe
        attempt<-attempt + 1
        BOLD_data_table <- data.frame(matrix(nrow=0, ncol=2))
      }
      #*****************************downloading from BOLD finished*******************************************************************************#

      if(nrow(BOLD_data_table)>0){

        summary_log<-paste0("Download Complete: cleaning BOLD data table - ", Sys.time(), " - ", genera_list[genera_list_loop_counter])
        (print(paste0("Download Complete: cleaning BOLD data table - ", genera_list[genera_list_loop_counter])))
        #Adding to the summary log for the running of the program and the results
        write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

        #Calling the clean function to clean the contents of the data_table
        BOLD_data_table_cleaned<-as.data.frame(data_clean(BOLD_data_table,genera_list[genera_list_loop_counter], seq_min, seq_max))

        #check if there is any data left after cleaning
        if(nrow(BOLD_data_table_cleaned) < 1){
          print(paste0("There is no useable data from BOLD download - " ,genera_list[genera_list_loop_counter]))
        }

      }else{

        print("BOLD Download Error: failed to populate data table or there is no data")
        #this ensures that any other data remaining inside this table is wiped
        BOLD_data_table_cleaned<-NULL
      }
    }#End of if BOLD_database
    #*****************************Begin NCBI section *******************************************************************************#

    #*****************************downloading from NCBI, attempt at least 3 times************************************************************************************************#

    if (NCBI_database==TRUE){

      summary_log<-paste0("Attempting NCBI Download - ", Sys.time(), " - ", genera_list[genera_list_loop_counter])
      (print(paste0("Attempting NCBI Download - ", genera_list[genera_list_loop_counter])))
      #Adding to the summary log for the running of the program and the results
      write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

      if (NCBI_flag == 1){

        #Create a subfolder to keep all of the initial fastas NCBI
        NCBI_folder_str<-paste0(main_file_folder, "/NCBI")
        dir.create(NCBI_folder_str)

        NCBI_flag = 0

      }

      #This is the section where I will build the search string for NCBI. If it is equal to the preset value then build it with the if, if it isn't then use it as the user inputted value
      if(is.null(search_str)){

        #This line builds the search string if no alternative was entered in as an attribute
        ncbi_search_str<-paste0("(", genera_list[genera_list_loop_counter],"[ORGN]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])")

      }else{

        ncbi_search_str<-search_str[genera_list_loop_counter]

      }

      #Using the below line to stop an error from occurring based on the NCBI system
      httr::set_config(httr::config(http_version = 1))

      #initialize the attempt loop counter, set the error occurrence flag to false and initialize an empty dataframe
      attempt<- 1
      NCBI_data_table <- data.frame(matrix(nrow=0, ncol=2))

      while(attempt < 4)
      {

        if(attempt == 2)
        {
          # This is setting the internet wait time to a higher value than the 60 sec default. We needed to make it at least 120 as we were having timeout errors.

          (print(paste0("In the attempt 2 of the error catching while loop with counter...", attempt)))
          Sys.sleep(0.1)
          options(timeout = 120)

        }
        else if (attempt == 3)
        {
          (print(paste0("In the attempt 3 of the error catching while loop with counter...", attempt)))
          Sys.sleep(0.1)
          options(timeout = 240)

        }

        NCBI_data_table<- tryCatch({
          #Call the NCBI_download function send it to the location for the output file
          result<-as.data.frame(NCBI_download(genera_list[genera_list_loop_counter], NCBI_folder_str, main_file_folder, ncbi_search_str))
        },
        error = function(e){
          print('Timeout Error, retrying with greater timeout limit')

        })

        #exit loop if error didn't occur
        if(typeof(NCBI_data_table) != 'character')
        {
          break
        }

        #increment the counter and reset the dataframe
        attempt<-attempt + 1
        NCBI_data_table <- data.frame(matrix(nrow=0, ncol=2))
      }


    #*****************************downloading from NCBI finished*******************************************************************************#

      if(nrow(NCBI_data_table)>0){
        summary_log<-paste0("Download Complete: cleaning NCBI data table - ", genera_list[genera_list_loop_counter])
        (print(paste0("Download Complete: cleaning NCBI data table - ", genera_list[genera_list_loop_counter])))
        #Adding to the summary log for the running of the program and the results
        write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

        #Calling the clean function to clean the contents of the data_table
        NCBI_data_table_cleaned<-as.data.frame(data_clean(NCBI_data_table,genera_list[genera_list_loop_counter], seq_min, seq_max))

        #check if theres any data left after cleaning
        if(nrow(NCBI_data_table_cleaned) < 1){
          print(paste0("There is no useable data from NCBI download - ",genera_list[genera_list_loop_counter]))
        }

      }
      else{
        print("NCBI Download Error: failed to populate data table or there is no data")
        #this ensures that any other data remaining inside this table is wiped
        NCBI_data_table_cleaned<-NULL
      }
    }# end of if NCBI_database

    if (is.null(BOLD_data_table_cleaned)){
      total_data_table<-NCBI_data_table_cleaned
    }else if (is.null(NCBI_data_table_cleaned)){
      total_data_table<-BOLD_data_table_cleaned
    }else if (nrow(BOLD_data_table_cleaned)==0){
      total_data_table<-NCBI_data_table_cleaned
    }else if (nrow(NCBI_data_table_cleaned)==0){
      total_data_table<-BOLD_data_table_cleaned
    }else if(nrow(BOLD_data_table_cleaned)>0 && nrow(NCBI_data_table_cleaned)>0){
      # Create the total data table for this taxa first we need to make the name values equal
      total_data_table<-rbind(BOLD_data_table_cleaned, NCBI_data_table_cleaned)
    }

    # This if is only adding this data to the total output file if there was data downloaded
    if(!is.na(nrow(total_data_table)) && nrow(total_data_table)>0){

      #Need to make a summary log file here.
      summary_log<-c(paste0("Number of records - ", nrow(total_data_table)))
      unique_sp_list<-paste(unique(total_data_table$Species),collapse=",")
      summary_log<-c(summary_log,paste0("Species - ",unique_sp_list))
      unique_markers_list<-paste(unique(total_data_table$Gene), collapse = ",")
      summary_log<-c(summary_log,paste0("Molecular markers - ", unique_markers_list))
      summary_log<-c(summary_log,paste0("Ending time - ", Sys.time(), " - ",genera_list[genera_list_loop_counter]))

      #Adding to the summary log for the running of the program and the results
      write.table(summary_log ,file=paste0(table_folder,"/A_Summary.txt"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\n", append=TRUE)

      #The next lines are making sure that the output is in a standard format for the other associated scripts
      #Removing the first row off the total_data_table
      total_data_table= total_data_table[-1,]

      #Initialize the output file for all of the molecular markers
      write.table(total_data_table ,file=paste0(table_folder,"/A_Total_Table.tsv"), na="", row.names=FALSE, col.names=FALSE, quote = FALSE,sep="\t", append=TRUE)

    }#Closing off the if statement checking to see if there was data for this taxa

  }# Closing off the main loop going through each of the taxa in the original list

} #Closing the main function for the BOLD/NCBI download
