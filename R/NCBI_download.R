#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'
#' @import rentrez

#********************************************Main program section***********************************************
############################################ NCBI FUNCTION #########################################################
NCBI_download <- function(target_genera, NCBI_folder_str, main_file_folder, search_str)
{

#******************************************** NCBI entrez DOWNLOAD ***********************************************

  (print(paste0("Here is the value of the search string...", search_str)))

  #completing the search using rentrez search function
  search_loop_target_taxa <- entrez_search(db="nuccore", term=search_str, use_history=TRUE)

  #building the output file name for the BOLD fasta
  NCBI_out_name<-paste0(NCBI_folder_str,"/",target_genera,"_NCBI_download.txt")

  #initializing the output variable from the loop
  NCBI_download_data<-NULL
  #This for loop is getting the records from the history on the NCBI folder from the 1 to the final (search_loop_target_taxa[2]) 10000 at a time
  for( seq_start in seq(0,as.numeric(search_loop_target_taxa[2]),1000)){
    print(paste0("Downloading 1000 ",target_genera , " NCBI records starting at...",seq_start, " of ", search_loop_target_taxa[2], " - ",Sys.time()))
    records <- entrez_fetch(db="nuccore", web_history=search_loop_target_taxa$web_history, rettype="gb", retmax=1000, retstart=seq_start)
    NCBI_download_data<-rbind(NCBI_download_data,records)
    #Writing fasta vector to a NCBI fasta file to read back in and then reformat the header and combine the BOLD and NCBI
    write(NCBI_download_data ,file=NCBI_out_name, append=FALSE)
  }


 #************************************************************ PRODUCE TABLE FROM NCBI DATA DOWNLOAD ***********************************

  # load in the data file in to the data for the lines
  NCBI_download_data_file<-readLines(NCBI_out_name)

  #Here the if statement checks to see if there is data in the file info vector. The first column represents the size fo the file of interest
  if(length(NCBI_download_data_file)!=0){

    # Get the line numbers with the end of each entry i the file
    record_ends<-grep(pattern = "//", NCBI_download_data_file)

    #initialize the starting position of the record to 1
    record_start=1

    #initialize the data frame to hold the total data_table
    NCBI_download_data_table<-data.frame(DB=character(), record_accession=character(), record_accession=character(), record_species=character(), DB=character(), record_gene=character(), record_sequence=character(), stringsAsFactors=FALSE)

    #Loop through each of the records in the main file to build a fasta formatted sequence
    for (NCBI_record_count in 1:length(record_ends)){

      #get all of the lines from the main file that belong to this record
      record<-NCBI_download_data_file[record_start:record_ends[NCBI_record_count]]

      record_accession<-record[grep(pattern = "ACCESSION", record)]

      #get the unique identifier in the accession number
      record_accession<-record[grep(pattern = "ACCESSION", record)]
      if (length(record_accession)!=0){
        record_accession<-gsub("ACCESSION","",record_accession)
        record_accession<-gsub(" ","",record_accession)
      }else{
        record_accession<-""
      }

      #get the species information with the ORGANISM identifier. The trims removes the leading white space
      record_species<-record[grep(pattern = "ORGANISM", record)]
      if (length(record_species)!=0){
        record_species<-gsub("ORGANISM","",record_species)
        record_species<-trimws(record_species)
      }else{
        record_species<-""
      }

      #Get the gene region. note some GenBank records have multiple instances of the gene identifier. I am only taking the first here.
      record_gene_p<-record[grep(pattern = "/product=", record)]
      record_gene_g<-record[grep(pattern = "/gene=", record)]
      if (length(record_gene_p)!=0){
        record_gene<-gsub("/product=","",record_gene_p[1])
        record_gene<-gsub(" ","",record_gene)
        record_gene<-gsub("\"","",record_gene)
      }else if (length(record_gene_g)!=0){
        record_gene<-gsub("/gene=","",record_gene_g[1])
        record_gene<-gsub(" ","",record_gene)
        record_gene<-gsub("\"","",record_gene)
      }else {
        record_gene<-""
      }

      #get all of the lines with sequence data
      record_sequence_start<-grep(pattern = "ORIGIN", record)
      
      #first if is ignoring any record that has many lines for the ORIGIN
      if(length(record_sequence_start)>10){
        record_sequence<-""
      } else if(length(record_sequence_start)!=0){
        #get all of the lines from the main file that belong to this record
        record_sequence<-record[record_sequence_start[1]:length(record)]
        record_sequence<-gsub("\\d","",record_sequence)
        record_sequence<-gsub(" ","",record_sequence)
        record_sequence<-paste0(record_sequence,collapse="")
        record_sequence<-gsub("//","",record_sequence)
        record_sequence<-gsub("ORIGIN","",record_sequence)
      } else {
        record_sequence<-""
      }

      #set the starting line for the next iteration of the loop to get the records
      record_start=record_ends[NCBI_record_count]+1

      #Here assemble tab delimited record for output to the data file
      data_record<-as.data.frame(cbind("GenBank", record_accession, record_accession, record_species, "GenBank", record_gene, record_sequence))

      NCBI_download_data_table<-rbind(NCBI_download_data_table, data_record)

    }#Closing the loop for each of the records in the main NCBI XML file

  } else{

    NCBI_download_data_table=NULL

  }#closing the if statement making sure that there is data in the downloaded file

  return(NCBI_download_data_table)

}
