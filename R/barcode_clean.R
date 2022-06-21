#Roxygen2 Documentation:

#' @export
#'
#' @title DNA Barcode Clean
#'
#' @author Robert G. Young and Jarrett D. Phillips
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
#' @param dist_model This is the genetic distance model of nucleotide substitution (as implemented through ape). The default model is "raw", which corresponds to simple p-distance. Other available models are JC69 (Jukes-Cantor, 1969), K80 (Kimura-2-Paramter), and F81 (Felenstein, 1981)
#' @param AGCT_only This indicates if records with characters other than AGCT are kept, the default is TRUE. TRUE removes records with non-AGCT FALSE is accepting all IUPAC characters
#' @param statistic This is the desired statistic to be resampled. By default, the DNA barcode gap (barcode_gap) is computed. Other statistics are the minimum interspecific distance (min_inter) and the maximum intraspecific distance (max_intra).
#' @param subsample_prop This is the subsample proportion used for bootstrapping, which should be between than 0 and 1 exclusive.
#' @param replicate_size This is is number of bootstrap replications. This value should be set to at least 1000. The default is 10000.
#' @param replacement This indicates sampling with replacement or sampling without replacement. The default is TRUE, indicating sampling with replacement.
#' @param conf_level This is the confidence level used for interval estimation. The default is 0.95, indicating 95% confidence.
#' @param data_folder This variable can be used to provide a location for the MSA fasta files to be cleaned. The default value is set to NULL where the program will prompt the user to select the folder through point-and-click.
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
#'

#********************************************Main program section***********************************************
##################################### Main FUNCTION ##############################################################
barcode_clean <- function(AA_code= c("invert", "vert", "std"),
                          dist_model = c("raw", "JC69", "K80", "F81"),
                          AGCT_only = TRUE,
                          statistic = c("barcode_gap", "min_inter", "max_intra"),
                          subsample_prop = NULL,
                          replicate_size = 10000,
                          replacement = TRUE,
                          conf_level = 0.95,
                          data_folder = NULL){

  #AA_code="invert"
  #AGCT_only = TRUE
  #data_folder = NULL

# Codes include 'std', 'vert', 'invert', 'NULL' skips the AA clean section
# AGCT_only TRUE is on and FALSE is accepting all IUPAC characters


  if (is.null(data_folder)){

    # prompting to choose the folder location of the working directory with the input file to run the program
    n <- substr(readline(prompt="Choose the folder location where your input files are located. Hit enter key to continue..."),1,1)
    #Get the directory
    Work_loc<-readpath()

  }else{

    Work_loc = data_folder

  }

#set the format for the date for all operating systems
Sys.setlocale("LC_TIME", "C")

# Current Date - for file naming use
date <- sub("-", "", sub("-", "", Sys.Date()))

#creating a log file
current_time<-as.character(Sys.time())
current_time<-gsub(" ","",current_time)
current_time<-gsub(":","",current_time)
log_file_name<- paste0("A_Clean_File_",current_time)
log_header<- paste0("DNA_Clean_Log_File - File name = ", log_file_name, " - AA code = ", AA_code, " - AGCT only = ", AGCT_only)


#Making the amino acid translation codes into numbers for the ape package.
if (AA_code == "vert"){
  AA_code = 2
}else if (AA_code == "invert"){
  AA_code = 5
}else if (AA_code == "std"){
  AA_code = 1
}else {
  (AA_code == 0)
}

#outputing the header to the log file
write.table(log_header, file=paste0(Work_loc,"/",log_file_name,".dat"),append=FALSE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

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

#************************** LOOPING THROUGH EACH OF THE FILES IN THE TARGET DIRECTORY *******************************************************************

for(h in 1:length(file_name)){

  file_sep<-paste0("********** ", file_name[h], " - ", Sys.time(), " **********")

  #outputting the file separator to the log file
  write.table(file_sep,file=paste0(Work_loc,"/",log_file_name,".dat"),append=TRUE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="")

  #initializing the variable to remove sequences for this loop
  seq_to_remove=NULL

  # load in the MSA
  Seq_file<-data.frame(readLines(path_list[h]))

  #using the read in file and transforming it to a two column table with header and sequence using the fasta_to_table function
  Seq_file_data_frame<-fasta_to_table(Seq_file)

  #creating a main separated file of the loaded in fasta to use throughout the script
  Seq_file_data_frame_total<-data.frame(do.call("rbind", strsplit(as.character(Seq_file_data_frame[,1]), "|", fixed = TRUE)))
  Seq_file_data_frame<-cbind(as.vector(Seq_file_data_frame[,1]), as.data.frame(Seq_file_data_frame_total), as.vector(Seq_file_data_frame[,2]))

  #Add column to the Seq_file_data_frame to hold the results of all of the checks
  Flags<-c("Flags")
  Seq_file_data_frame[,Flags]<-"-"

  #Give headers for the new data frame
  colnames(Seq_file_data_frame)<-c("Header", "ID", "Accession", "Genus", "Species", "BIN_OTU", "Gene", "Sequence", "Flags")

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

    #Getting rid of sequences with non AGCT characters
    no_AGCT_seq<-subset(Seq_file_data_frame,grepl("[^AGCTagct-]",Seq_file_data_frame$Sequence))

    #Here I get entries which don't have a flag
    flag_subset<-subset(Seq_file_data_frame,Seq_file_data_frame$Flags == "-")

    #checking to see if there are records without a flag
    if((as.numeric(nrow(flag_subset))>0) && (as.numeric(nrow(no_AGCT_seq))>0) ){

      #Getting the records with - in Flag and have characters other than AGCT and are listed in no_AGCT_seq
      flag_subset<-flag_subset[flag_subset$Header %in% no_AGCT_seq$Header,]

      #Adding the results of AA check to the Seq_file_data_frame Flags column
      Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "non_AGCT"

    }

    #taking the row names of the reduced data set
    seq_to_remove<-unique(c(seq_to_remove,no_AGCT_seq$Header))

    #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
    no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]

    if(nrow(no_outliers_dataset)>0){

      #Then determine the number of each taxa present in the data set by using a freq table
      ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

      #Rename the last column to AGCT
      colnames(ft_sp)<-c("Species","AGCT")

      #Here I am merging the log_df with the results after the AGCT clean
      log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

    }
  }

  #**************************** Removing sequences with stop codons *********************************

  if(AA_code!=0 && nrow(Seq_file_data_frame)>2){

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
      flag_subset<-flag_subset[flag_subset$Header %in% row.names(AA_string),]

      #Adding the results of AA check to the Seq_file_data_frame Flags column
      Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Stop_Codon"

    }

    #taking the row names of the identified records with stop codons and adding it to the seq_to_remove list
    seq_to_remove<-c(seq_to_remove,rownames(AA_string))

    #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
    no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]

    if(nrow(no_outliers_dataset)>0){

      #Then determine the number of each taxa present in the data set by using a freq table
      ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

      #Rename the last column to AA
      colnames(ft_sp)<-c("Species","AA")

      #Here I am merging the log_df with the results after the AA clean
      log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)
    }
  }


  #**************************** Outlier based cleaning for Genus *********************************

  if(nrow(Seq_file_data_frame)>2){

    #using the ape function to obtain the distance matrix
    dist_matrix<-suppressWarnings(dist.dna(Seq_file_DNAbin, model = dist_model, variance = FALSE, gamma = FALSE, pairwise.deletion = TRUE, base.freq = NULL, as.matrix = TRUE))

    #Making the contents of the matrix numeric
    dist_matrix<-apply(dist_matrix, 2, as.numeric)

    #This is adding the row names back to the matrix as the above line removed it.
    row.names(dist_matrix)<-colnames(dist_matrix)

    #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
    no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),]
    no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove)]

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
      flag_subset<-flag_subset[flag_subset$Header %in% seq_to_remove_outlier,]

      #Adding the results of AA check to the Seq_file_data_frame Flags column
      Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Genus_Outlier"

    }

    #Add the seq_to_remove_outlier to the total seq_to_remove list
    seq_to_remove<-c(seq_to_remove,seq_to_remove_outlier)

    #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
    no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]

    if(nrow(no_outliers_dataset)>0){

      #Then determine the number of each taxa present in the data set by using a freq table
      ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

      #Rename the last column to Outlier
      colnames(ft_sp)<-c("Species","Genus_Outlier")

      #Here I am merging the log_df with the results after the Outlier clean
      log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

    }

  }#Closing the if to check if there are more than two records to evaluate so the outlier calcs can work

  #**************************** Species outlier to report in the records output only - NOT REDUCING THE DATASET JUST REPORTING *********************************

  if(nrow(Seq_file_data_frame)>2){

    # Get a list of the species in the genus
    Species<-unique(no_outliers_dataset$Species)

    #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
    no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),]
    no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove)]

    #loop through the species
    for (species_outlier_list_counter in 1:length(Species)){

      #Get the records for the outlier loop species - I can use the no_outliers_dataset
      #because I still have it from the Genus outlier check
      outlier_loop_species_records<-no_outliers_dataset[no_outliers_dataset$Species==Species[species_outlier_list_counter],]

      if(nrow(outlier_loop_species_records)>2) {

        #Get the rows of the no outgroup distance matrix for the species of interest
        outlier_loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% outlier_loop_species_records$Header),]

        #Get the columns with the species of interest from the rows of the species of interest.
        outlier_loop_species_dist_matrix_within <- outlier_loop_species_dist_matrix[,(colnames(outlier_loop_species_dist_matrix) %in% outlier_loop_species_records$Header)]

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
          flag_subset<-flag_subset[flag_subset$Header %in% seq_to_remove_outlier,]

          #Adding the results of check to the Seq_file_data_frame Flags column
          Seq_file_data_frame$Flags[Seq_file_data_frame$Header %in% c(flag_subset$Header)]<- "Species_Outlier"

        }

        #Add the seq_to_remove_outlier to the total seq_to_remove list
        seq_to_remove<-c(seq_to_remove,seq_to_remove_outlier)

      }# End of the loop to check to see if there are records for the species in question

    }# end of loop going through species in the genus

      #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
      no_outliers_dataset<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]

    if(nrow(no_outliers_dataset)>0){

      #Then determine the number of each taxa present in the data set by using a freq table
      ft_sp<-as.data.frame(table(no_outliers_dataset$Species))

      #Rename the last column to Outlier
      colnames(ft_sp)<-c("Species","Species_Outlier")

      #Here I am merging the log_df with the results after the Outlier clean
      log_df<-merge(log_df,ft_sp,by=1, all.x=TRUE)

    }

  }#end of the species more than 2 records check

  #**************************** Barcode Gap analysis section *********************************

  if(nrow(Seq_file_data_frame)>2){

    #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
    no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),]
    no_outliers_dist_matrix <- no_outliers_dist_matrix[, !(colnames(no_outliers_dist_matrix) %in% seq_to_remove)]

    #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
    barcode_gap_data_frame<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]

    #Establishing a list of the genera present in the data set, this should only be one
    Genera<-unique(barcode_gap_data_frame$Genus)

    #Add barcode gap reporting columns to the species output reporting table log_df
    barcode_gap_columns<-c("Intraspecific", "Interspecific", "Barcode_Gap", "Barcode_Gap_Value", "Estimate_Bias", "Estimate_SE",  "Estimate_CI_Lower", "Estimate_CI_Upper")
    log_df[,barcode_gap_columns]<-"-"

    #Initializing the storage data frame for the barcode gap results
    barcode_gap_results<-NULL

    #Initializing the plotting data frame
    plot_out<-data.frame()

    #Initializing the plotting flag so that the first loop is equal to the first colum of records
    plot_flag=0

    #This section will subset the data in the loop into each of the Genera present in the data set one at a time
    for (genera_list_counter in 1:length(Genera)){

      #Get the subset of the data with the Genera for this loop
      loop_genus_records<-barcode_gap_data_frame[barcode_gap_data_frame$Genus==Genera[genera_list_counter],]

      if(nrow(loop_genus_records)>1){

        #Sort the data frame by Species
        loop_genus_records<- loop_genus_records[order(loop_genus_records$Species),]

        #Get a unique species list
        Species<-unique(loop_genus_records$Species)

          for (species_list_counter in 1:length(Species)){

            #Get the records for the loop species
            loop_species_records<-loop_genus_records[loop_genus_records$Species==Species[species_list_counter],]

            #for species with more than one records so that the within is able to be calculated

            if (nrow(loop_species_records)>1 && length(Species)>1){

              #Get the rows of the target species from the dist matrix and then get the columns from the selected columns
              loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% loop_species_records$Header),]
              loop_species_dist_matrix_within <- loop_species_dist_matrix[,(colnames(loop_species_dist_matrix) %in% loop_species_records$Header)]

              #Now get comparisons between the loop species and all other records
              loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% loop_species_records$Header),]
              loop_species_dist_matrix_between <- loop_species_dist_matrix[,!(colnames(loop_species_dist_matrix) %in% loop_species_records$Header)]

              ##### Resampling to calculate barcode gap standard error (SE) #####
              # perform resampling - Added by Jarrett

              # select desired statistic
              statistic <- match.arg(statistic)

              # preallocate vector of resamples
              boot_samples <- numeric(replicate_size)

              # resample subsample_size genetic distances with or without replacement replicate_size times
              for (i in 1:replicate_size) {
                if (replacement == TRUE) { # bootstrapping
                  intra_boot <- sample(loop_species_dist_matrix_within, size = ceiling(subsample_prop * length(loop_species_dist_matrix_within)), replace = TRUE)
                  inter_boot <- sample(loop_species_dist_matrix_between, size = ceiling(subsample_prop * length(loop_species_dist_matrix_between)), replace = TRUE)
                } else { # subsampling
                  intra_boot <- sample(loop_species_dist_matrix_within, size = ceiling(subsample_prop * length(loop_species_dist_matrix_within)), replace = FALSE)
                  inter_boot <- sample(loop_species_dist_matrix_between, size = ceiling(subsample_prop * length(loop_species_dist_matrix_between)), replace = FALSE)
                }

                if (statistic == "barcode_gap") {
                  # bootstrapped barcode gap
                  boot_samples[i] <- min(inter_boot) - max(intra_boot)
                  # observed sample barcode gap
                  stat_obs <- min(loop_species_dist_matrix_between) - max(loop_species_dist_matrix_within)
                } else if (statistic == "min_inter") {
                  # bootstrapped minimum intraspecific distance
                  boot_samples[i] <- min(inter_boot)
                  # observed sample minimum interspecfic distance
                  stat_obs <- min(loop_species_dist_matrix_between)
                } else { # max_intra
                  # bootstrapped minimum intraspecific distance
                  boot_samples[i] <- max(intra_boot)
                  # observed sample maximum intraspecific distance
                  stat_obs <- max(loop_species_dist_matrix_within)
                }
              }

              # calculate  bootstrap mean
              stat_boot_mean <- mean(boot_samples)

              # calculate bootstrap standard error
              stat_boot_se <- sd(boot_samples)

              # calculate percentile CI
              stat_boot_ci <- quantile(boot_samples, c((1 - conf_level) / 2, (1 + conf_level) / 2))

              #Getting the maximum within species distance
              loop_species_dist_matrix_within<-max(loop_species_dist_matrix_within)

              #Getting the minimum between distance
              loop_species_dist_matrix_between<-min(loop_species_dist_matrix_between)

              #calculating the gap
              if(loop_species_dist_matrix_within<loop_species_dist_matrix_between){

                loop_species_barcode_gap<-"YES"

              }else{

                loop_species_barcode_gap<-"NO"

              }

            } else{

              loop_species_dist_matrix_within<-"NA"
              loop_species_dist_matrix_between<-"NA"
              loop_species_barcode_gap<-"NA"

            }#closing the if else checking if there is more than one species record

          #add the results of the species barcode gap within or intraspecific to the log_df
          log_df$Intraspecific[log_df$Species %in% Species[species_list_counter] ]<-loop_species_dist_matrix_within

          #add the results of the species barcode gap between or interspecific to the log_df
          log_df$Interspecific[log_df$Species %in% Species[species_list_counter] ]<-loop_species_dist_matrix_between

          #add the results of the species barcode gap check to the log_df
          log_df$Barcode_Gap[log_df$Species %in% Species[species_list_counter] ]<-loop_species_barcode_gap

          #add the results of the species barcode gap calculation to the log_df
          log_df$Barcode_Gap_Value[log_df$Species %in% Species[species_list_counter] ] <- loop_species_dist_matrix_between - loop_species_dist_matrix_within

          #add results of the bootstrap SE
          log_df$Estimate_SE[log_df$Species %in% Species[species_list_counter] ] <- stat_boot_se

          # calculate bootstrap bias
          stat_boot_bias <- stat_boot_mean - stat_obs

          #add results of the bootstrap SE
          log_df$Estimate_Bias[log_df$Species %in% Species[species_list_counter] ] <- stat_boot_bias

          #add results of the lower bootstrap CI endpoint
          log_df$Estimate_CI_Lower[log_df$Species %in% Species[species_list_counter] ] <- stat_boot_ci[1]

          #add results of the upper bootstrap CI endpoint
          log_df$Estimate_CI_Upper[log_df$Species %in% Species[species_list_counter] ] <- stat_boot_ci[2]

          # plot sampling distribution
          par(mfrow = c(1, 2))

          hist(boot_samples) # histogram
          abline(v = stat_boot_mean, lty = 2)
          qqnorm(boot_samples) # QQ plot
          qqline(boot_samples)

        }#closing the loop through the unique species in the genus

      }#closing if checking if there is more than one record in the genus

    }#end of the more than 2 records check

  }#closing loop through unique genera

  #***************************************** OUTPUT ********************************************

  #Getting the data ready to be written to file

  #Printing the distance matrix to a file
  write.table(no_outliers_dist_matrix,file=paste0(Work_loc,"/",file_name[h],"_dist_matrix.dat"),append=FALSE,na="",row.names = TRUE, col.names=TRUE, quote = FALSE,sep="\t")

  #Writing the data to the data table
  write.table(Seq_file_data_frame,file=paste0(Work_loc,"/",file_name[h],"_data_table.dat"),append=FALSE,na="",row.names = FALSE, col.names=TRUE, quote = FALSE,sep="\t")

  if(!is.null(seq_to_remove)){

    #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
    no_outliers_dataset<-no_outliers_dataset[!(no_outliers_dataset$Header %in% seq_to_remove),]

    #removing the columns with data other than the header and the sequence data to construct the fasta file.
    no_outliers_dataset<-cbind(no_outliers_dataset$Header, no_outliers_dataset$Sequence)

    #outputting the fasta reduced file
    write.table(no_outliers_dataset,file=paste0(Work_loc,"/",file_name[h],"_no_outliers.fas"),append=FALSE,na="",row.names = FALSE, col.names=FALSE, quote = FALSE,sep="\n")

    #Adding the results to the log file
    #outputting the file separator to the log file
    suppressWarnings(write.table(log_df,file=paste0(Work_loc,"/",log_file_name,".dat"),append=TRUE,row.names = FALSE, col.names=TRUE, quote = FALSE,sep="\t"))

  }

  print(log_df)

} #end of file loop

print(paste("Start time... ",start_time," and end time... ",Sys.time()))

} #End of the function

