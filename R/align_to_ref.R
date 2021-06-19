######################### MAIN FUNCTION FOR ALIGNMENT WITH MAFFT #####################
#Roxygen2 Documentation:

#' @export
#'
#' @title Alignment and trimming of a multiple sequence alignment to a reference sequence
#'
#' @author Robert G. Young
#'
#' @description
#' This function takes a FASTA file with target sequences and aligns them against a reference sequence submitted to the program.
#' The output is an aligned fasta file that is trimmed to the length of the reference sequence. Sequences without full coverage (records having sequences with leading or trailing gaps) are removed.
#' Records with characters other than IUPAC are also removed. Finally, internal gaps are removed from the sequence based on the submitted multiple sequence alignment percent
#' coverage of the character position as provided in the pigl argument supplied by the user.
#'
#' @details
#' User Input:
#'     1. A file folder location with the fasta files that need to be aligned and trimmed using the supplied reference sequence. Please note that any and all fasta files (named *.fas) in this folder will be analyzed.
#'     2. A reference sequence file with a sequence or MSA with all sequences having the same length.
#'     3. The location of the MAFFT executable file (https://mafft.cbrc.jp/alignment/software/)
#'
#' Some example formats for the running of the function are...
#' align_to_ref(pigl=0.75)
#' align_to_ref(pigl=0.95, op=10)
#' align_to_ref(pigl=0)
#'
#' @param pigl This is the percent internal gap loop argument. This provides a percent that will remove records causing internal gaps if more than the percent value assigned
#' to this argument is reached. If this value is set to 0 then internal gaps are not removed. The default for this value is 0.95.
#'
#' @param op This is the gap opening penalty for the use of MAFFT. The higher the value the larger penalty in the alignment. The default for this value is set to 1.53 which is the
#' default value in the MAFFT program. For alignment of highly conserved regions where no gaps are expected this should be set to a much higher number and 10 is recommended for coding regions like the COI-5P.
#'
#' @returns Output:
#' 1. In the submitted file folder location there will be a log file titled “MAFFT_log”.
#' 2. The sequence output files from this script are placed into two subfolders. These folders are in the submitted file location where the fasta files of interest are located.
#'    The two folders created are MAFFT and MAFFT_trimmed. In the MAFFT folder there will be files with name of the files in the submitted file folder appended with “_MAFFT”.
#'    The MAFFT_trimmed file will contain files with the same naming convention as the files in the submitted folder and appended with “MAFFT_trimmed”.
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young, R. G., Gill, R., Gillis, D., Hanner, R. H. (Submitted June 2021). Molecular Acquisition, Cleaning, and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal.
#'
#' @seealso
#' auto_seq_download()
#' create_fastas()
#' barcode_clean()
#'
#' @import utils
#'


# This function aligns a sequence file to a reference sequence and then trims to the beginning and end of the reference sequence
##################################### align_to_ref FUNCTION ##############################################################
align_to_ref <- function(pigl=0.95, op=1.53){

  #pigl (per internal gap in loop) is the fractional value of the percent of the sequence that needs to have gap characters at a particular position before that gap is removed.

  start_wd<-getwd()

  # prompting to choose the file folder
  n <- substr(readline(prompt="Choose the folder location where your fasta files to be aligned are located. Hit enter key to continue..."),1,1)
  #This line is requesting the file folder location for the fasta files of interest
  Work_loc<-readpath()

  # Choose the fasta file with the reference sequence to align the new sequence
  n <- substr(readline(prompt="Choose your fasta reference file (note this must be a trimmed file with all sequences of the same length and no leading or trailing gap characters). Hit enter key to continue..."),1,1)
  reference_MSA<-file.choose()

  #Make a file folder to contain all results from the program
  main_file_folder<-paste0(Work_loc, "/MAFFT")
  dir.create(main_file_folder)
  main_file_folder<-paste0(Work_loc, "/MAFFT_trimmed")
  dir.create(main_file_folder)

  #Make the log file for the running of the script
  write.table("",file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

  #Get all files located in this folder.
  #Initiating the list of files in the folders
  path_list <- list.files(path=Work_loc, pattern = "*[.][Ff][Aa][Ss]$", full.names = TRUE)
  file_list <- list.files(path=Work_loc, pattern = "*[.][Ff][Aa][Ss]$")
  file_names <- sub("\\..*","",as.vector(file_list))

  # prompting to choose the folder location for the MAFFT .bat program
  n <- substr(readline(prompt="Choose the folder location where the MAFFT (.bat file) is located. Hit enter key to continue..."),1,1)
  MAFFT_location<-readpath()
  setwd(MAFFT_location)

  #initializing the final total fasta output variable
  output_total_fasta<-NULL

  # loop through all of the files in the file folder. b
  for(i in 1:length(file_list)){

    #printing the file name and time stamp to console
    name_and_time_stamp<-paste0(path_list[i], " at ", Sys.time())
    print(name_and_time_stamp)

    #Writing the file name and date and time to the log file
    write.table(name_and_time_stamp,file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #This was the area where I create the location and name of the output fasta file
    output_MSA<-paste0(Work_loc,"/MAFFT/",sub("\\..*","_MAFFT.fas",as.vector(file_list[i])))

    #run the MAFFT alignment program
    cmd_mafft_string<- paste("mafft --add", as.vector(path_list[i]), "--op", op, reference_MSA, ">", output_MSA)

    #Writing the command string to the console and log file
    write.table(paste0("MAFFT command... ", cmd_mafft_string),file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #running the MAFFT program by submitting the built string above to a command shell window
    MAFFT_output<-system(cmd_mafft_string, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL, show.output.on.console = FALSE, minimized = FALSE, invisible = TRUE, timeout = 0)

    #output the MAFFT specifics to the log file
    write.table(MAFFT_output, file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #Read in the new alignment file and trim to target region, also read in the reference sequence alignment and remove these form the final alignment

    # load in the reference MSA and the output MSA data file in
    reference_fasta_file<-data.frame(readLines(reference_MSA))
    output_fasta_file<-data.frame(readLines(output_MSA))

    #use the fasta_to_table function and take the read in fasta file and place it in a table for use
    reference_fasta_table<-fasta_to_table(reference_fasta_file)
    output_fasta_table<-fasta_to_table(output_fasta_file)

    #making the output fasta table sequence column upper case
    output_fasta_table[,2]<-toupper(output_fasta_table[,2])

    #using the MSA_to_matrix function to take the two column (Header, sequence) table and make it a matrix with each nucleotide position in a column
    output_fasta_matrix<-as.data.frame(MSA_data_matrix(output_fasta_table))

    #printing specifics to the console and log file
    term_log_info<-paste0("Length of the fasta file reported in nucleotides ", nchar(as.character(output_fasta_matrix[1,2])), " and the number of records ", nrow(output_fasta_matrix))
    print(term_log_info)

    #output the MAFFT specifics to the log file
    write.table(term_log_info, file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #identify the start and stop positions of the sequences from the reference file in the total alignment
    #first extract the columns from the ouput_fasta_matrix which were in the reference fasta file

    #get the reference file headers
    seq_in_reference<-reference_fasta_table[,1]

    #subset the output_fasta
    output_fasta_matrix_reference_subset<-subset(output_fasta_matrix,output_fasta_matrix[,1]==seq_in_reference)

    #getting the trim locations for the output file
    forward_reverse_pos<-forward_reverse_trim(output_fasta_matrix_reference_subset)

    #removing the reference sequences from the output fasta matrix
    output_fasta_matrix<-as.data.frame(subset(output_fasta_matrix,output_fasta_matrix[,1]!=seq_in_reference))

    #deleting the forward and end gaps on either side of the reference sequence
    if(forward_reverse_pos[2]!=ncol(output_fasta_matrix)&forward_reverse_pos[1]!=1){
      output_fasta_matrix<-as.data.frame(output_fasta_matrix[,-c(3:(forward_reverse_pos[1]-1),(forward_reverse_pos[2]+1):ncol(output_fasta_matrix))])
    }else if(forward_reverse_pos[1]==1){
      output_fasta_matrix<-as.data.frame(output_fasta_matrix[,-c((forward_reverse_pos[2]+1):ncol(output_fasta_matrix))])
    }else if(forward_reverse_pos[2]==ncol(output_fasta_matrix)){
      output_fasta_matrix<-as.data.frame(output_fasta_matrix[,-c(3:(forward_reverse_pos[1]-1))])
    }

    #Collapsing all of the columns into a single column
    output_fasta_matrix[,2]<-as.data.frame(apply(output_fasta_matrix[,3:ncol(output_fasta_matrix)],1,paste,collapse=""))

    # Remove records with no sequence remaining after trimming (so sequences with all - and/or N)
    output_fasta_matrix<-as.data.frame(subset(output_fasta_matrix,nchar(gsub("[-N]","",output_fasta_matrix[,2]))!=0))

    # If the internal gap is set to 0 then skip the removal of internal gaps for non coding genes
    if(pigl !=0 ){
      # Deleting internal gaps where there is greater than the pigl percentage of gaps
      output_fasta_matrix<-as.data.frame(remove_internal_gaps(output_fasta_matrix, pigl))
    }

    #Collapsing all of the columns in to a single column and building the output two column table for output
    output_fasta_matrix[,2]<-as.data.frame(apply(output_fasta_matrix[,3:ncol(output_fasta_matrix)],1,paste,collapse=""))
    output_fasta_matrix<-as.data.frame(output_fasta_matrix[,-c(3:ncol(output_fasta_matrix))])

    #subsetting the data.frame to remove the sequences with N's
    output_fasta_matrix<-as.data.frame(subset(output_fasta_matrix,(nchar(gsub("[^N]+","",output_fasta_matrix[,2])))==0))

    #Get the number of characters in the alignment
    characters_in_alignment<-nchar(as.character(output_fasta_matrix[1,2]))

    #Remove leading and trailing gaps from sequences in the alignment
    output_fasta_matrix[,2]<-as.data.frame(gsub("^-+|-+$", "", output_fasta_matrix[,2]))

    #Removing sequences with fewer than the target number of nucleotides in the final alignment
    output_fasta_matrix<-as.data.frame(subset(output_fasta_matrix, nchar(as.character(output_fasta_matrix[,2]))==characters_in_alignment))

    # If the internal gap is set to 0 then skip the removal of internal gaps for non coding genes
    if(pigl !=0 ){

      #subsetting the data.frame to remove the sequences with gaps
      output_fasta_matrix<-as.data.frame(subset(output_fasta_matrix,(nchar(gsub("[^-]+","",output_fasta_matrix[,2])))==0))

    }

    #printing specifics to the console and log file
    term_log_info<-paste0("Length of the aligned and trimmed multiple sequence alignment in nucleotides ", nchar(as.character(output_fasta_matrix[1,2])), " and the number of records ", nrow(output_fasta_matrix))
    print(term_log_info)

    #output the MAFFT specifics to the log file
    write.table(term_log_info, file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #output the A division line in the log file
    write.table("********************************************************************************", file=paste0(Work_loc,"/MAFFT_log.txt"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #order the sequences by genus_species name
    #also copying the header and pulling out the Genus Header
    gen_header<-gsub("^.*?\\|","",output_fasta_matrix[,1])
    gen_header<-gsub("^.*?\\|","",gen_header)
    gen_header<-gsub("^.*?\\|","",gen_header)
    gen_header<-gsub("\\|.*", "", gen_header)

    #also copying the header and pulling out the Species Header
    sp_header<-gsub("^.*?\\|","",output_fasta_matrix[,1])
    sp_header<-gsub("^.*?\\|","",sp_header)
    sp_header<-gsub("^.*?\\|","",sp_header)
    sp_header<-gsub("^.*?\\|","",sp_header)
    sp_header<-gsub("\\|.*", "", sp_header)
    gen_sp_header<-paste(gen_header,sp_header,sep="_")

    #Adding the genus_species column and then ordering the martix by genus_species and then removing the column
    output_fasta_matrix<-as.data.frame(cbind(output_fasta_matrix,gen_sp_header))
    output_fasta_matrix<-as.data.frame(output_fasta_matrix[order(output_fasta_matrix[,3]),])
    output_fasta_matrix<-as.data.frame(output_fasta_matrix[,-3])

    #This was the area where I create the location and name of the output fasta file
    output_MSA<-paste0(Work_loc,"/MAFFT_trimmed/",sub("\\..*","_MAFFT_Trimmed.fas",as.vector(file_list[i])))

    #outputting the MAFFT trimmed file
    write.table(output_fasta_matrix,file=output_MSA,append=FALSE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    print("Output is located in the target directory in the subfolders MAFFT and MAFFT_trimmed")

  }
  setwd(start_wd)

}

