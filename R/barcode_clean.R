#Roxygen2 Documentation:

#' @export
#'
#' @title DNA Barcode Clean
#'
#' @author Robert G. Young
#'
#' @description
#' Takes an input fasta file and identifies genus level outliers and species outliers based on the 1.5 x greater than the interquartile range.
#' It also, if selected, checks the sequence using amino acid translation and has the option to eliminate sequences that have non-IUPAC codes.
#' Finally, the program calculates the barcode gap for the species in the submitted dataset.
#'
#' @details
#' Input: A file folder with one or more fasta files of interest
#'
#' @examples
#' \dontrun{
#' barcode_clean(),
#' barcode_clean(AA_code = "vert", AGCT_only = TRUE),
#' barcode_clean(AA_code = "vert")
#' }
#'
#' @param AA_code This is the amino acid translation matrix (as implemented through ape) used to check the sequences for stop codons. The following codes are available std, vert, invert, F. The default is invert.
#' @param AGCT_only This indicates if records with characters other than AGCT are kept, the default is TRUE. TRUE removes records with non-AGCT FALSE is accepting all IUPAC characters
#' @param data_folder This variable can be used to provide a location for the MSA fasta files to be cleaned. The default value is set to NULL where the program will prompt the user to select the folder through point-and-click.
#' @param gen_outliers This is the variable to indicate if the user would like to remove suspected genus sequence record outliers using 1.5X the genetic distance. If set to TRUE genus level outliers will be removed. If FALSE this will not occur. Default TRUE.
#' @param sp_outliers This is the variable to indicate if the user would like to remove suspected species sequence record outliers using 1.5X the genetic distance. If set to TRUE species level outliers will be removed. If FALSE this will not occur. Default TRUE.
#' @param dist_model This is the model of nucleotide evolution that the ape program will use (see ape documentation for options. Default is "raw"
#' @param replicates This is the number of replicates that the bootstrapping will perform. Note: more replicates will take longer. Default is 1000
#' @param replacement This indicates that the replacement of MSA nucleotide columns will be replaced in the random resampling. Default is set to TRUE
#' @param conf_level This is a percentage of the initial MSA nucleotide length. When set to 1 the bootstrapped resampling will have the same length as the initial MSA. Default is set to 1
#' @param numCores This is the number of cores that the user would like to use where multithreading is available. Default is set to 1, indicating only a single thread will be used.
#'
#' @returns
#' Output:
#' A single log file for the running of the function with the name A_Clean_File_YYYY-DD-TTTTTTTT.
#' The function will also output three files for each fasta file submitted. The first is the distance matrix that was calculated and used to assess the DNA barcode gaps.
#' This file is named the same as the input file with dist_table.dat appended to the end of the name. The second file is the total data table file which provides a table
#' of all submitted records for each data set accompanied with the results from each section of the analysis. This file is named the same as the input fasta with data_table.dat appended
#' to the end, Finally, a fasta file with all outliers and flagged records removed is generated for each input fasta file. This output file is named the same as the input fasta with no_outlier.fas appended to the end.
#'
#' @references
#' <https://github.com/rgyoung6/MACER>
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'
#' @seealso
#' auto_seq_download()
#' create_fastas()
#' align_to_ref()
#'
#' @import ape
#' @import stats
#' @import utils
#' @import ggplot2
#' @import parallel
#' @import pbapply
#' @import grDevices
#' @import png

#********************************************Main program section***********************************************
##################################### Main FUNCTION ##############################################################
barcode_clean <- function(AA_code="invert", AGCT_only = TRUE, data_folder = NULL, gen_outliers = TRUE, sp_outliers = TRUE, dist_model = "raw", replicates = 1000, replacement = TRUE,conf_level = 1,  numCores = 1){

#  AA_code=0
#  AGCT_only = TRUE
#  data_folder = NULL
#  gen_outliers = FALSE
#  sp_outliers = FALSE
#  dist_model = "raw"
#  replicates = 1000
#  replacement = TRUE
#  conf_level = 1
#  numCores = 1

#  library(ape)
#  library(stats)
#  library(utils)
#  library(ggplot2)
#  library(parallel)
#  library(pbapply)
#  library(grDevices)
#  library(png)


  #Get the initial working directory
  start_wd <- getwd()
  on.exit(setwd(start_wd))

  dateStamp <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M"), "_")

  if (is.null(data_folder)){

    # prompting to choose the folder location of the working directory with the input file to run the program
    n <- substr(readline(prompt="Choose a file in the folder location where your input files are located. Hit enter key to continue..."),1,1)
    #Get the directory
    Work_loc<-dirname(file.choose())

  }else{

    Work_loc = data_folder

  }

  setwd(Work_loc)

  #set the format for the date for all operating systems
  Sys.setlocale("LC_TIME", "C")

  # Current Date - for file naming use
  date <- sub("-", "", sub("-", "", Sys.Date()))

  #creating a log file
  current_time<-as.character(Sys.time())
  current_time<-gsub(" ","",current_time)
  current_time<-gsub(":","",current_time)
  log_file_name<- paste0("A_Clean_File_",current_time)
  log_header<- paste0("DNA_Clean_Log_File - File name = ", log_file_name, " - AA code = ", AA_code, " - AGCT only = ", AGCT_only,
                      " dist_model = ", dist_model, " replicates = ", replicates, " replacement = ", replacement," conf_level = ", conf_level, " numCores = ", numCores)

  #Making the amino acid translation codes into numbers for the ape package.
  if (is.null(AA_code)){
    AA_code = 0
  }else if (AA_code == "vert"){
    AA_code = 2
  }else if (AA_code == "invert"){
    AA_code = 5
  }else if (AA_code == "std"){
    AA_code = 1
  }else if (AA_code == "moldPlus"){
    AA_code = 4
  }else if (AA_code == "yeast"){
    AA_code = 3
  }else {
    AA_code = 0
  }

  start_time=Sys.time()
  print(paste("Start time... ",start_time))

  #Get all files located in this folder.
  #Initiating the list of files in the folders
  path_list <- list.files(path=Work_loc, pattern = "*[.][Ff][Aa][Ss]$", full.names = TRUE)
  file_list <- list.files(path=Work_loc, pattern = "*[.][Ff][Aa][Ss]$")
  file_name <- sub("\\..*","",as.vector(file_list))

  #Initialize the no_outliers_dist_matrix
  no_outliers_dist_matrix=""
  no_outliers_dataset=""

  #outputing the header to the log file
  write.table(log_header, file=paste0(Work_loc,"/",log_file_name,".dat"),append=FALSE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

  #************************** LOOPING THROUGH EACH OF THE FILES IN THE TARGET DIRECTORY *******************************************************************

  for(h in 1:length(file_name)){

    file_sep<-paste0("********** ", file_name[h], " - ", Sys.time(), " **********")

    print(paste0(file_sep))

    #initializing the variable to remove sequences for this loop
    seq_to_remove=NULL

    #initialize the variable to print to the log file for each analysed file
    loop_flag = 1

    # load in the MSA
    Seq_file<-data.frame(readLines(path_list[h]))

    #using the read in file and transforming it to a two column table with header and sequence using the fasta_to_table function
    Main_Seq_file_data_frame<-fasta_to_table(Seq_file)

    #creating a main separated file of the loaded in fasta to use throughout the script
    Main_Seq_file_data_frame_total<-data.frame(do.call("rbind", strsplit(as.character(Main_Seq_file_data_frame[,1]), "|", fixed = TRUE)))
    Main_Seq_file_data_frame<-cbind(as.vector(Main_Seq_file_data_frame[,1]), as.data.frame(Main_Seq_file_data_frame_total), as.vector(Main_Seq_file_data_frame[,2]))

    #Add column to the Main_Seq_file_data_frame to hold the results of all of the checks
    Flags<-c("Flags")
    Main_Seq_file_data_frame[,Flags]<-"-"

    #Give headers for the new data frame
    colnames(Main_Seq_file_data_frame)<-c("Header", "ID", "Accession", "Genus", "Species", "BIN_OTU", "Gene", "Sequence", "Flags")

    #Remove Duplicates
    Main_Seq_file_data_frame<- Main_Seq_file_data_frame[!duplicated(Main_Seq_file_data_frame$Header),, drop=FALSE]

    #For each genus in the data file
    unique_genera_list <- unique(Main_Seq_file_data_frame$Genus)

    for(unique_genera in 1:length(unique_genera_list)){

      print(paste("Genus ", unique_genera_list[unique_genera], " - ",unique_genera, " of ", length(unique_genera_list)))
      #Get the dataframe with the records from the loop genus
      Seq_file_data_frame <- Main_Seq_file_data_frame[Main_Seq_file_data_frame$Genus == unique_genera_list[unique_genera],, drop=FALSE]

      if(nrow(Seq_file_data_frame)>1){

        # freq table of species
        ft_sp<-as.data.frame(table(Seq_file_data_frame$Species))

        #Add in the species list and initial number of sequences to populate the log_df
        log_df<- data.frame(Species=as.character(ft_sp[,1]),Initial=as.integer(ft_sp[,2]))

        #Get the list of records without an accession number indicating that they are only in BOLD to save and add back after dereplicating the
        #records based on accession that are present in both GenBank and BOLD
        no_accession<-as.data.frame(subset(Seq_file_data_frame, Seq_file_data_frame$Accession=="NA"))

        #Get all records with an accession number
        yes_accession<-as.data.frame(subset(Seq_file_data_frame, Seq_file_data_frame$Accession!="NA"))

        #dereplicate the yes_accession data set based on accession column
        remove_duplicates<-as.data.frame(subset(yes_accession, !duplicated(yes_accession$Accession)))

        #combine the no accession and dereplicated yes accession data sets
        Seq_file_data_frame<-rbind(remove_duplicates, no_accession)

        # freq table of species to use to add to log file with number of sequences after dereplication
        ft_sp<-as.data.frame(table(Seq_file_data_frame$Species))

        #Rename the last column to Dereplicate
        colnames(ft_sp)<-c("Species","Dereplicate")

        #Here I am merging the log_df with the results after the remove duplicates
        log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

        #Getting the sequence data and converting it to DNAbin
        #The t is transposing the matrix. strsplit the second column containing the sequence data in to columns for each character
        suppressWarnings(Seq_file_DNAbin<-do.call("rbind", strsplit(as.character(Seq_file_data_frame$Sequence), "", fixed = TRUE)))

        # Making the row names of the new matrix equal to the headers for the fasta file
        rownames(Seq_file_DNAbin)<-Seq_file_data_frame$Header

        # Making the file of class DNAbin for use with ape functions
        Seq_file_DNAbin<-as.DNAbin(Seq_file_DNAbin)

        #**************************** Removing sequences with non AGCT characters *********************************

        if(AGCT_only==TRUE && nrow(Seq_file_data_frame)>2){

          print(paste0("Removing records with non-AGCT characters at ", Sys.time()))

          #Getting rid of sequences with non AGCT characters
          no_AGCT_seq<-subset(Seq_file_data_frame,grepl("[^AGCTagct-]",Seq_file_data_frame$Sequence))

          #Here I get entries which don't have a flag
          flag_subset<-subset(Seq_file_data_frame,Seq_file_data_frame$Flags == "-")

          #checking to see if there are records without a flag
          if((as.numeric(nrow(flag_subset))>0) && (as.numeric(nrow(no_AGCT_seq))>0) ){

            #Getting the records with - in Flag and have characters other than AGCT and are listed in no_AGCT_seq
            flag_subset<-flag_subset[flag_subset$Header %in% no_AGCT_seq$Header,, drop=FALSE]

            #Adding the results of AA check to the Seq_file_data_frame Flags column
            Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "non_AGCT"

          }

          #taking the row names of the reduced data set
          seq_to_remove<-unique(c(seq_to_remove,no_AGCT_seq$Header))

          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),, drop=FALSE]

          if(nrow(no_outliers_dataset)>0){

            #Then determine the number of each taxa present in the data set by using a freq table
            ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

            #Rename the last column to AGCT
            colnames(ft_sp)<-c("Species","AGCT")

            #Here I am merging the log_df with the results after the AGCT clean
            log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

          }
        }else{
          #Add a column with NA for the AGCT spot
          log_df$AGCT <- "NA"
        }

        #**************************** Removing sequences with stop codons *********************************

        if(AA_code!=0 && nrow(Seq_file_data_frame)>2){

          print(paste0("Removing recrods with stop codons at ", Sys.time()))

          # Using ape function take Seq_file_DNAbin and translate into AA
          # Code 1 is standard code, 2 is vertebrate mitochondrial, 5 is invert mitochondrial and is an argument at the beginning of the
          AA <- trans(Seq_file_DNAbin, code = AA_code, codonstart = 0)

          # Take the AAbin object format and make it into a matrix
          AA_matrix<-as.matrix(as.character(AA))

          #Collapsing all AA into a single string
          AA_string<-as.data.frame(apply(AA_matrix,1,paste,collapse=""))

          #Getting the records with stop codons
          AA_string<-subset(AA_string,grepl("\\*",AA_string[,1]))

          #Here I get entries which don't have a flag from the main Seq_file_data_frame
          flag_subset<-subset(Seq_file_data_frame,Seq_file_data_frame$Flags == "-")

          #checking to see if there are records without a flag
          if((as.numeric(nrow(flag_subset))>0) && (as.numeric(nrow(AA_string))>0) ){

            #Getting the records with - in Flag and have a stop codon
            flag_subset<-flag_subset[flag_subset$Header %in% row.names(AA_string),, drop=FALSE]

            #Adding the results of AA check to the Seq_file_data_frame Flags column
            Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Stop_Codon"

          }

          #taking the row names of the identified records with stop codons and adding it to the seq_to_remove list
          seq_to_remove<-c(seq_to_remove,rownames(AA_string))

          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),, drop=FALSE]

          if(nrow(no_outliers_dataset)>0){

            #Then determine the number of each taxa present in the data set by using a freq table
            ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

            #Rename the last column to AA
            colnames(ft_sp)<-c("Species","AA")

            #Here I am merging the log_df with the results after the AA clean
            log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)
          }
        }else{
          #Add a column with NA for the AGCT spot
          log_df$AA <- "NA"
        }

        #*************** Build the distance matrix for the cleaned sequence data *************************

        if(nrow(Seq_file_data_frame)>2){

          print(paste0("Building the distance matrix for the molecular data set at ", Sys.time()))

          #using the ape function to obtain the distance matrix
          if (dist_model == "raw" || dist_model == "JC69" || dist_model == "K80" || dist_model == "F81") {
            dist_matrix<-suppressWarnings(dist.dna(Seq_file_DNAbin, model = dist_model, variance = FALSE, gamma = FALSE, pairwise.deletion = TRUE, base.freq = NULL, as.matrix = TRUE))
          } else {
            stop(paste("No calculation of the distance matrix as the model selected (",dist_model, ") is not an accepted model. Please select an appropriate model and start the function again"))
          }

          #Making the contents of the matrix numeric
          dist_matrix<-apply(dist_matrix, 2, as.numeric)

          #Making the diagonals of the symetric matrix NA so that there are no self comparisons
          diag(dist_matrix)<-NA

          #This is adding the row names back to the matrix as the above line removed it.
          row.names(dist_matrix)<-colnames(dist_matrix)

          #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
          no_outliers_dist_matrix <- dist_matrix[!(rownames(dist_matrix) %in% seq_to_remove),, drop=FALSE]
          no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove), drop=FALSE]

        }#output is the no_outliers_dist_matrix

        #**************************** Outlier based cleaning for Genus *********************************

        if(gen_outliers == TRUE && nrow(Seq_file_data_frame)>2){

          print(paste0("Genus outlier assessment at ", Sys.time()))

          #Getting the inter quartiles
          dist_matrix_quant <- quantile(no_outliers_dist_matrix, na.rm=TRUE)

          #Get the inter quartile range which can be used to multiple 1.5 times. This added to the upper quartile will give the outlier bounds
          dist_matrix_inter_quant = as.numeric((dist_matrix_quant[4])+(1.5*((as.numeric(dist_matrix_quant[4])-as.numeric(dist_matrix_quant[2])))))

          #initiate the seq_to_remove_outlier variable
          seq_to_remove_outlier<-NULL

          for (p in 1:nrow(no_outliers_dist_matrix)){
            if (length(na.omit(no_outliers_dist_matrix[,p]))>1){

              if((mean(as.numeric(na.omit(no_outliers_dist_matrix[,p])))>=dist_matrix_inter_quant) && (median(as.numeric(na.omit(no_outliers_dist_matrix[,p])))>=dist_matrix_inter_quant) && (dist_matrix_inter_quant!=0)){
                seq_to_remove_outlier<-c(seq_to_remove_outlier,rownames(no_outliers_dist_matrix)[p])
              }
            }
          }

          #Here I get entries which don't have a flag
          flag_subset<-subset(Seq_file_data_frame,Seq_file_data_frame$Flags == "-")

          #checking to see if there are records without a flag
          if((as.numeric(nrow(flag_subset))>0) && (as.numeric(length(seq_to_remove_outlier))>0) ){

            #Getting the records with - in Flag and are statistical outliers by genetic distance
            flag_subset<-flag_subset[flag_subset$Header %in% seq_to_remove_outlier,, drop=FALSE]

            #Adding the results of AA check to the Seq_file_data_frame Flags column
            Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Genus_Outlier"

          }

          #Add the seq_to_remove_outlier to the total seq_to_remove list
          seq_to_remove<-c(seq_to_remove,seq_to_remove_outlier)

          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),, drop=FALSE]

          if(nrow(no_outliers_dataset)>0){

            #Then determine the number of each taxa present in the data set by using a freq table
            ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

            #Rename the last column to Outlier
            colnames(ft_sp)<-c("Species","Genus_Outlier")

            #Here I am merging the log_df with the results after the Outlier clean
            log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

          }

        }else{#Closing the if to check if there are more than two records to evaluate so the outlier calcs can work

          #Add a column with NA for the Genus_Outlier spot
          log_df$Genus_Outlier <- "NA"
        }

        #**************************** Species outlier to report in the records output only - NOT REDUCING THE DATASET JUST REPORTING *********************************

        if(sp_outliers == TRUE && nrow(Seq_file_data_frame)>2){

          print(paste0("Species outlier assessment at ", Sys.time()))

          # Get a list of the species in the genus
          Species<-unique(no_outliers_dataset$Species)

          #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
          no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),, drop=FALSE]
          no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove), drop=FALSE]

          #loop through the species
          for (species_outlier_list_counter in 1:length(Species)){

            #Get the records for the outlier loop species - I can use the no_outliers_dataset
            #because I still have it from the Genus outlier check
            outlier_loop_species_records<-no_outliers_dataset[no_outliers_dataset$Species==Species[species_outlier_list_counter],, drop=FALSE]

            if(nrow(outlier_loop_species_records)>2) {

              #Get the rows of the no outgroup distance matrix for the species of interest
              outlier_loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% outlier_loop_species_records$Header),, drop=FALSE]

              #Get the columns with the species of interest from the rows of the species of interest.
              outlier_loop_species_dist_matrix_within <- outlier_loop_species_dist_matrix[,(colnames(outlier_loop_species_dist_matrix) %in% outlier_loop_species_records$Header), drop=FALSE]

              #Getting the inter quartiles for only the species to species comparisons
              dist_matrix_quant <- quantile( outlier_loop_species_dist_matrix_within , na.rm=TRUE)

              #Get the inter quartile range which can be used to multiple 1.5 times. This added to the upper quartile will give the outlier bounds
              dist_matrix_inter_quant = as.numeric((dist_matrix_quant[4])+(1.5*((as.numeric(dist_matrix_quant[4])-as.numeric(dist_matrix_quant[2])))))

              #initiate the seq_to_remove_outlier variable
              seq_to_remove_outlier<-NULL

              #getting the records with the genetic distance greater than the calculated
              for (p in 1:nrow(outlier_loop_species_dist_matrix_within)){
                if (length(na.omit(outlier_loop_species_dist_matrix_within[,p]))>1){

                  if((mean(as.numeric(na.omit(outlier_loop_species_dist_matrix_within[,p])))>=dist_matrix_inter_quant) && (median(as.numeric(na.omit(outlier_loop_species_dist_matrix_within[,p])))>=dist_matrix_inter_quant) && (dist_matrix_inter_quant!=0)){
                    seq_to_remove_outlier<-c(seq_to_remove_outlier,rownames(outlier_loop_species_dist_matrix_within)[p])
                  }
                }
              }

              #Here I get entries which don't have a flag
              flag_subset<-subset(Seq_file_data_frame,Seq_file_data_frame$Flags == "-")

              #checking to see if there are records without a flag
              if((as.numeric(nrow(flag_subset))>0) && (as.numeric(length(seq_to_remove_outlier))>0) ){

                #Getting the records with - in Flag and are statistical outliers by genetic distance
                flag_subset<-flag_subset[flag_subset$Header %in% seq_to_remove_outlier,, drop=FALSE]

                #Adding the results of check to the Seq_file_data_frame Flags column
                Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Species_Outlier"

              }

              #Add the seq_to_remove_outlier to the total seq_to_remove list
              seq_to_remove<-c(seq_to_remove,seq_to_remove_outlier)

            }# End of the loop to check to see if there are records for the species in question

          }# end of loop going through species in the genus

          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),, drop=FALSE]

          if(nrow(no_outliers_dataset)>0){

            #Then determine the number of each taxa present in the data set by using a freq table
            ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

            #Rename the last column to Outlier
            colnames(ft_sp)<-c("Species","Species_Outlier")

            #Here I am merging the log_df with the results after the Outlier clean
            log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

          }

        }else{#end of the species more than 2 records check

          #Add a column with NA for the Species_Outlier spot
          log_df$Species_Outlier <- "NA"
        }
        #**************************** Setting up the Barcode Gap analysis section *********************************

        if(nrow(Seq_file_data_frame)>2){

          #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
          no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),]
          no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove), drop=FALSE]

          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          barcode_gap_data_frame<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),, drop=FALSE]

          #Sort the data frame by Species
          barcode_gap_data_frame<- barcode_gap_data_frame[order(barcode_gap_data_frame$Species),, drop=FALSE]

          #Get a unique species list
          Species<-unique(barcode_gap_data_frame$Species)

          #Add barcode gap reporting columns to the species output reporting table log_df
          barcode_gap_columns<-c("Target_Records", "Target_Haplotypes", "Genus_Records", "Genus_Haplotypes", "Num_Overlap_Haploptype", "Overlap_Haploptype","Intraspecific", "Interspecific", "Barcode_Gap", "Num_Barcode_Gap_Overlap_Records", "Barcode_Gap_Overlap_Records", "Num_Barcode_Gap_Overlap_Taxa", "Barcode_Gap_Overlap_Taxa", "Barcode_Gap_Result", "Bootstrap")
          log_df[,barcode_gap_columns] <- "-"

        }

        #**************************** Setting up the bootstrapping matricies *********************************

        if(nrow(Seq_file_data_frame)>2){

          #If the replicates is greater than 0 then need to bootstrap
          if(replicates != 0){

            print(paste0("Bootstrapping section at ", Sys.time()))

            #Initialize the list
            rep_matricies <- vector("list", replicates)

            #Reduce the target genus data set to only include unique taxa and sequence records for use in the bootstrapping analysis
            bootstrap_records <- barcode_gap_data_frame[!duplicated(barcode_gap_data_frame[,c("Genus", "Species", "Sequence")]),c("Header", "Sequence")]

            #Break the fasta file into unique columns for the sequences
            bootstrap_records<-as.data.frame(cbind(bootstrap_records$Header, t(as.data.frame(strsplit(as.character(bootstrap_records$Sequence), "")))))

            #Sort the dataframe by headers
            bootstrap_records<-bootstrap_records[order(bootstrap_records[,1]),, drop=FALSE]

            #Make the Headers column the rownames
            row.names(bootstrap_records)<-bootstrap_records[,1]
            bootstrap_records <- bootstrap_records[,-1]
            colnames(bootstrap_records)<-c(1:ncol(bootstrap_records))

            #Make bootstrap_records a data frame
            bootstrap_records<-as.data.frame(bootstrap_records)

            #Get a list of species for this genus
            list_of_species<-unique(Seq_file_data_frame$Species)
            boot_species_results<-as.data.frame(cbind(list_of_species,boot_results=0))
            boot_species_results_final<-as.data.frame(list_of_species)
            row.names(boot_species_results_final)<-boot_species_results_final[,1]
            boot_species_results_final<-boot_species_results_final[,-1, drop=FALSE]

            bootstrap_barcode_gap <- function(i){

              #Sample the target species records to create a resample matrix with replacement
              boot_loop_records <- sample(colnames(bootstrap_records), size = ceiling(ncol(bootstrap_records)*conf_level), replace = replacement)

              #Create the within matrix using the resample with replacement
              boot_loop_matrix_temp <- as.data.frame(bootstrap_records[,match(boot_loop_records, colnames(bootstrap_records))], check.names = FALSE, drop = FALSE)

              # Making the file of class DNAbin for use with ape functions
              boot_loop_matrix_temp_DNAbin<-as.DNAbin(as.matrix(boot_loop_matrix_temp))

              #using the ape function to obtain the distance matrix
              if (dist_model == "raw" || dist_model == "JC69" || dist_model == "K80" || dist_model == "F81") {
                boot_loop_dist_matrix<-as.data.frame(suppressWarnings(dist.dna(boot_loop_matrix_temp_DNAbin, model = dist_model, variance = FALSE, gamma = FALSE, pairwise.deletion = TRUE, base.freq = NULL, as.matrix = TRUE)))
              } else {
                stop(paste("No calculation of the distance matrix as the model selected (",dist_model, ") is not an accepted model. Please select an appropriate model and start the function again"))
              }

              #For each species calculate the presence of the barcode gap and then save to a dataframe and return the dataframe
              for (species_count in 1:length(list_of_species)){
                between_values<-boot_loop_dist_matrix[grepl( list_of_species[species_count] , rownames( boot_loop_dist_matrix ) ) , !grepl( list_of_species[species_count] , colnames( boot_loop_dist_matrix ) ) ,drop=FALSE]
                within_values<-boot_loop_dist_matrix[grepl( list_of_species[species_count] , rownames( boot_loop_dist_matrix ) ) , grepl( list_of_species[species_count] , colnames( boot_loop_dist_matrix ) ) ,drop=FALSE]
                if(!is.null(between_values) & !is.null(within_values)){
                  if(ncol(between_values)>0 & nrow(within_values)>0 ){
                    if((as.numeric(min(between_values))- as.numeric(max(within_values)))>0){
                      #There is a barcode gap so add one to the results
                      boot_species_results[boot_species_results[,1]==list_of_species[species_count],2]<-as.numeric(boot_species_results[boot_species_results[,1]==list_of_species[species_count],2]) + 1
                    }
                  }
                }
              }

              return(boot_species_results[,2])

            } #End of bootstrap function

            if (numCores == 1){
              boot_species_results_final <- cbind(boot_species_results_final, pblapply(seq_len(replicates), bootstrap_barcode_gap))
            }else if(numCores >1){
              boot_species_results_final <- cbind(boot_species_results_final, pblapply(seq_len(replicates), bootstrap_barcode_gap, cl = numCores))
            }
          }
        }

        #*********************** Begin looping through the species in the genus ***************************

        if(nrow(Seq_file_data_frame)>2){

          for (species_list_counter in 1:length(Species)){

            #Get the records for the loop species
            loop_species_records<-barcode_gap_data_frame[barcode_gap_data_frame$Species==Species[species_list_counter],, drop= FALSE]

            #Initialize the variable to hold the bootstrap results
            loop_species_bootstrap = 0

            #for species with more than one records so that the within is able to be calculated
            if (nrow(loop_species_records)>1 && length(Species)>1){

              #Get the rows of the target species from the dist matrix and then get the columns from the selected columns
              loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% loop_species_records$Header),, drop=FALSE]
              loop_species_dist_matrix_within <- loop_species_dist_matrix[,(colnames(loop_species_dist_matrix) %in% loop_species_records$Header), drop=FALSE]

              #Now get comparisons between the loop species and all other records
              loop_species_dist_matrix_between <- loop_species_dist_matrix[,!(colnames(loop_species_dist_matrix) %in% loop_species_records$Header), drop=FALSE]

              #Get the number of records for the target species
              loop_species_target<-nrow(loop_species_dist_matrix_within)

              #Get the number of non target records in the genus
              loop_species_nontarget<-ncol(loop_species_dist_matrix_between)

              #Getting the maximum within species distance
              loop_species_dist_matrix_within<-round(max(loop_species_dist_matrix_within, na.rm=TRUE), digits = 4)

              #Getting the minimum between distance
              loop_species_dist_matrix_between_dist<-round(min(loop_species_dist_matrix_between, na.rm=TRUE), digits = 4)

              #calculating the gap
              loop_species_barcode_gap<-round(as.numeric(loop_species_dist_matrix_between_dist) - as.numeric(loop_species_dist_matrix_within), digits = 4)

              #Reporting the barcode gap overlapping taxa
              if(nrow(loop_species_dist_matrix_between)>1){
                #Get the min value of each column and place into a matrix row name and min value
                loop_species_barcode_gap_overlap_records<- as.data.frame(apply(loop_species_dist_matrix_between,2,min))
              }

              #Drop min values greater than the max between to get the overlapping records
              loop_species_barcode_gap_overlap_records <- loop_species_barcode_gap_overlap_records[loop_species_barcode_gap_overlap_records[,1, drop = FALSE] <= loop_species_dist_matrix_within,, drop=FALSE]

              #Get the overlapping unique haplotypes
              if (nrow(loop_species_barcode_gap_overlap_records) == 0){
                loop_species_overlap_haplotypes = "-"
              }else{
                loop_species_overlap_haplotypes <- paste(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% row.names(loop_species_barcode_gap_overlap_records), "Sequence"]), collapse = ", ")
              }

              if(nrow(loop_species_barcode_gap_overlap_records) != 0){ #Overlap records

                #Get the number of records overlapping
                loop_species_num_barcode_gap_overlap_records <- nrow(loop_species_barcode_gap_overlap_records)

                #Get the list of taxa from the records by first splitting the headers
                loop_species_num_barcode_gap_overlap_taxa <- data.frame(do.call("rbind", strsplit(as.character(row.names(loop_species_barcode_gap_overlap_records)), "|", fixed = TRUE)))[,c(3,4)]

                #Report the number of overlapping haplotypes
                loop_species_num_overlap_haplotypes <-length(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% row.names(loop_species_barcode_gap_overlap_records), "Sequence"]))

                #Check to see if the number of unique species overlap haplotypes is greater than 20
                if(length(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% row.names(loop_species_barcode_gap_overlap_records), "Sequence"]))<20){

                  loop_species_overlap_haplotypes <- paste(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% row.names(loop_species_barcode_gap_overlap_records), "Sequence"]), collapse = ", ")

                }else{

                  loop_species_overlap_haplotypes <- "Greater than 20"

                }

                #If to place a greater than 20 in the reporting variable if too many records
                if(nrow(loop_species_barcode_gap_overlap_records) < 20){

                  #Using the record names get all records.
                  loop_species_barcode_gap_overlap_records <- paste0(row.names(loop_species_barcode_gap_overlap_records), collapse = ", ")

                }else{

                  loop_species_barcode_gap_overlap_records <- "Greater than 20"

                }

                #Keeping the genus and species and creating a binomial name
                loop_species_num_barcode_gap_overlap_taxa <- paste0(loop_species_num_barcode_gap_overlap_taxa[,1]," ", loop_species_num_barcode_gap_overlap_taxa[,2])

                #Get the taxa
                loop_species_barcode_gap_overlap_taxa <- unique(loop_species_num_barcode_gap_overlap_taxa)

                #If to place a greater than 20 in the reporting variable if too many records

                #Using the record names get all records.
                loop_species_barcode_gap_overlap_taxa <- paste0(loop_species_barcode_gap_overlap_taxa, collapse = ", ")

                #Get the number of taxa
                loop_species_num_barcode_gap_overlap_taxa <- length(unique(loop_species_num_barcode_gap_overlap_taxa))

                #Collapsing all unique taxa into a single variable.
                loop_species_num_barcode_gap_overlap_taxa <- paste0(as.vector(unique(loop_species_num_barcode_gap_overlap_taxa)),collapse = ", ")

              }else{

                loop_species_num_barcode_gap_overlap_records = "-"
                loop_species_num_barcode_gap_overlap_taxa = "-"
                loop_species_barcode_gap_overlap_taxa = "-"
                loop_species_num_overlap_haplotypes <- "-"
                loop_species_overlap_haplotypes <- "-"
              }

              if(loop_species_barcode_gap <= 0){
                loop_species_result <- "NO"
              }else{
                loop_species_result <- "YES"
              }

              #Get the unique target sequences for the haplotype reporting
              loop_target_haplotypes <- length(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% loop_species_records$Header, "Sequence"]))

              #Get the unique off target genera sequences for the haplotype reporting
              loop_genus_haplotypes <- length(unique(barcode_gap_data_frame[!(barcode_gap_data_frame$Header) %in% loop_species_records$Header, "Sequence"]))

              #There was clearly more than one record as we made it here past the more than one species record if so that means that the two records or more were
              #dereplicated. So, then that would mean that the between species distance will be 0 and the only record would have this as the comparison
              #Create the variable to report
              if(replicates !=0){
                loop_species_bootstrap_result <- sum(as.numeric(as.character(boot_species_results_final[row.names(boot_species_results_final)==Species[species_list_counter],])))
              }else{
                loop_species_bootstrap_result = "-"
              }

            }else{ #Closing the if more than one species if

              loop_species_target <- "-"
              loop_target_haplotypes <- "-"
              loop_species_nontarget <- "-"
              loop_genus_haplotypes <- "-"
              loop_species_num_overlap_haplotypes <- "-"
              loop_species_overlap_haplotypes <- "-"
              loop_species_dist_matrix_within <- "-"
              loop_species_dist_matrix_between_dist <- "-"
              loop_species_barcode_gap <- "-"
              loop_species_num_barcode_gap_overlap_records <- "-"
              loop_species_barcode_gap_overlap_records <- "-"
              loop_species_num_barcode_gap_overlap_taxa <- "-"
              loop_species_barcode_gap_overlap_taxa <- "-"
              loop_species_result <- "-"
              loop_species_bootstrap_result <- "-"

            }#closing the if else checking if there is more than one species record

            #add the results of the species barcode gap within or intraspecific to the log_df
            log_df$Target_Records[log_df$Species %in% Species[species_list_counter] ]<-loop_species_target

            #Add the results of the target species haplotypes to the log df
            log_df$Target_Haplotypes[log_df$Species %in% Species[species_list_counter] ]<-loop_target_haplotypes

            #add the results of the species barcode gap within or intraspecific to the log_df
            log_df$Genus_Records[log_df$Species %in% Species[species_list_counter] ]<-loop_species_nontarget

            #Add the results of the genus haplotypes to the log df
            log_df$Genus_Haplotypes[log_df$Species %in% Species[species_list_counter] ]<-loop_genus_haplotypes

            #Add the number of overlapping exact sequence haplotypes
            log_df$Num_Overlap_Haploptype[log_df$Species %in% Species[species_list_counter] ] <- loop_species_num_overlap_haplotypes

            #Add the results of the target sequences that have the exact same sequence in the non target species in the genus
            log_df$Overlap_Haploptype[log_df$Species %in% Species[species_list_counter] ] <- loop_species_overlap_haplotypes

            #add the results of the species barcode gap within or intraspecific to the log_df
            log_df$Intraspecific[log_df$Species %in% Species[species_list_counter] ]<-loop_species_dist_matrix_within

            #add the results of the species barcode gap between or interspecific to the log_df
            log_df$Interspecific[log_df$Species %in% Species[species_list_counter] ]<-loop_species_dist_matrix_between_dist

            #add the results of the species barcode gap check to the log_df
            log_df$Barcode_Gap[log_df$Species %in% Species[species_list_counter] ]<-loop_species_barcode_gap

            #add the number of results of the species barcode gap overlapping records to the log_df
            log_df$Num_Barcode_Gap_Overlap_Records[log_df$Species %in% Species[species_list_counter] ]<-loop_species_num_barcode_gap_overlap_records

            #add the results of the species barcode gap overlapping records to the log_df
            log_df$Barcode_Gap_Overlap_Records[log_df$Species %in% Species[species_list_counter] ]<-loop_species_barcode_gap_overlap_records

            #add the results of the species barcode gap overlapping taxa to the log_df
            log_df$Num_Barcode_Gap_Overlap_Taxa[log_df$Species %in% Species[species_list_counter] ]<-loop_species_num_barcode_gap_overlap_taxa

            #add the results of the species barcode gap overlapping taxa to the log_df
            log_df$Barcode_Gap_Overlap_Taxa[log_df$Species %in% Species[species_list_counter] ] <- loop_species_barcode_gap_overlap_taxa

            #add the results of the species barcode gap check to the log_df
            log_df$Barcode_Gap_Result[log_df$Species %in% Species[species_list_counter] ]<-loop_species_result

            #add the results of the species bootstrapping
            log_df$Bootstrap[log_df$Species %in% Species[species_list_counter] ]<-loop_species_bootstrap_result

            #Get the row for this loop to output to the file and
            #Add the genus to the front of the species name being outputted
            to_print<-log_df[log_df$Species == Species[species_list_counter] ,]
            to_print$Species <-paste0(unique_genera_list[unique_genera], " ", to_print$Species)

            #Ensure that the values are characters so it will output properly
            to_print<-data.frame(lapply(to_print, as.character), stringsAsFactors=FALSE)

            #Replace all NA with -
            to_print[to_print=="NA"] <- "-"
            to_print[is.na(to_print)] <- "-"
            to_print[to_print=="numeric(0)"] <- "-"

            #only print to file for the first element in the loop
            if(loop_flag == 1){

              #outputting the file separator to the log file
              write.table(file_sep,file=paste0(Work_loc,"/",log_file_name,".dat"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="")
              suppressWarnings(write.table(to_print,file=paste0(Work_loc,"/",log_file_name,".dat"),append=TRUE,row.names = FALSE, col.names=TRUE, quote = FALSE,sep="\t"))

              loop_flag=0

            }else{

              #Adding the results to the log file
              suppressWarnings(write.table(to_print,file=paste0(Work_loc,"/",log_file_name,".dat"),append=TRUE,row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\t"))

            }

          }#closing the loop through the unique species in the genus

          #***************************************** OUTPUT ********************************************

          #Getting the data ready to be written to file

          #Printing the distance matrix to a file
          write.table(no_outliers_dist_matrix,file=paste0(Work_loc,"/",file_name[h],"_",unique_genera_list[unique_genera],"_dist_matrix.dat"),append=FALSE,na="",row.names = TRUE, col.names=NA, quote = FALSE,sep="\t")

          #Writing the data to the data table
          write.table(Seq_file_data_frame,file=paste0(Work_loc,"/",file_name[h],"_",unique_genera_list[unique_genera],"_data_table.dat"),append=FALSE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\t")

          if(!is.null(seq_to_remove)){

            #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
            no_outliers_dataset<-no_outliers_dataset[!(no_outliers_dataset$Header %in% seq_to_remove),]

            #removing the columns with data other than the header and the sequence data to construct the fasta file.
            no_outliers_dataset<-cbind(no_outliers_dataset$Header, no_outliers_dataset$Sequence)

            #outputting the fasta reduced file
            write.table(no_outliers_dataset,file=paste0(Work_loc,"/",file_name[h],"_no_outliers.fas"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

          }

        } #End of species in genus section

      }#End of the if more than two records in the genus loop dataset

    }#End of the genus loop

  } #end of file loop

  print(paste("Start time... ",start_time," and end time... ",Sys.time()))

} #End of the function
