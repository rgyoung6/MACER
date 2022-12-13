#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young RG, Gill R, Gillis D, Hanner RH (2021) Molecular Acquisition, Cleaning and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal 9: e71378. <https://doi.org/10.3897/BDJ.9.e71378>
#'
############################################ DATA TABLE CLEAN FUNCTION ####################################################
data_clean <- function(clean_data_table, target_genus, seq_min, seq_max)
{

  #Add column to the clean_data_table to hold the results of all of the checks
  Flags<-c("Flags")
  clean_data_table[,Flags]<-"-"

  #Adding an index to the beginning of the imported file to use as a indicator of what record I am at
  clean_data_table<-cbind(a=c(1:nrow(clean_data_table)),clean_data_table)

  #Name all of the columns to use as reference below because the names differ between the BOLD and NCBI downloads
  colnames(clean_data_table) <- c("uniqueID","DB", "ID", "Accession", "Genus_Species", "BIN_OTU", "Gene", "Sequence", "Flags")

  #Here I flag entries which don't have species
  flag_subset<-subset(clean_data_table$uniqueID,clean_data_table$Genus_Species == "")

  #Adding the results of the above check to the clean_data_table Flags column
  clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "No_Taxa"

  #Here I flag entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  #checking to see if there are records without a flag
  if(nrow(flag_subset>0)){

    #Flagging records with numbers in the species names
    flag_subset<-subset(flag_subset$uniqueID, grepl("\\d+",flag_subset$Genus_Species))

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Taxa_w_digits"
  }

  #Here I flag entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  #checking to see if there are records without a flag
  if(nrow(flag_subset>0)){

    #Flagging records with species names having punctuation
    flag_subset<-subset(flag_subset$uniqueID, grepl("[[:punct:]]+",flag_subset$Genus_Species))

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Taxa_w_punct"

  }

  #Function to count the number of spaces in the genus species names
  countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

  #Here I flag entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  #checking to see if there are records without a flag
  if(nrow(flag_subset>0)){

    #Flagging records with genus species names with more than one space indicating that there are more than two names
    flag_subset<-subset(flag_subset$uniqueID, as.numeric(countSpaces(as.character(flag_subset$Genus_Species)))>1)

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Taxa_name_issue"

  }

  #Getting full records from the main clean_data_table without flags
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  #checking to see if there are records without a flag and taking the table and splitting and adding genus and species to the end of the table
  if(nrow(flag_subset>0)){

    #Create a data frame with genus and species in different columns for the whole table subset with out flags
    gen_sp<- t(as.data.frame(strsplit(as.character(flag_subset$Genus_Species), " ",fixed=TRUE)))

    #Renaming the flag_subset column headers so I can rbind them later with the same column headers as the flagged records
    colnames(gen_sp)<-c("Genus","Species")

    #This is creating a data frame with the whole table subset without flags and the now separated two columns for the genus and species
    flag_subset<-cbind(flag_subset,gen_sp)

  }

  #Getting full records from the main clean_data_table with flags to then rbind the unflagged records
  flag_subset_flagged<-subset(clean_data_table,clean_data_table$Flags != "-")

  if(nrow(flag_subset_flagged)>0){

    #Adding two columns of Genus and Species at then and of the flag_subset_flagged records so that the table is equal to that of the unflagged records with the two new columns of Genus and Species
    gen_sp_headers<-c("Genus", "Species")
    flag_subset_flagged[,gen_sp_headers]<-"-"

    #taking the flagged and the unflagged records and combining them through rbind
    clean_data_table<-rbind(flag_subset,flag_subset_flagged)

  } else{

    clean_data_table<-flag_subset

  }

  #Reorder the data frame and drop the genus_species column
  clean_data_table<-as.data.frame(cbind(clean_data_table["uniqueID"], clean_data_table["DB"], clean_data_table["ID"], clean_data_table["Accession"], clean_data_table["Genus_Species"], clean_data_table["Genus"], clean_data_table["Species"], clean_data_table["BIN_OTU"], clean_data_table["Gene"], clean_data_table["Sequence"], clean_data_table["Flags"]))

  #Here I get entries which don't have a flag
  flag_subset<-as.data.frame(subset(clean_data_table, clean_data_table$Flags == "-"), drop=FALSE)

  if(nrow(flag_subset)>0){

    #Here I flag entries which don't have sequences
    flag_subset<-subset(flag_subset$uniqueID,flag_subset$Sequence=="")

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "No_Seq"

  }

  #Here I get entries which don't have a flag
  flag_subset<-as.data.frame(subset(clean_data_table, clean_data_table$Flags == "-"), drop=FALSE)

  if(nrow(flag_subset)>0){

    #Here I flag entries which have sequences either below the desired size
    flag_subset<-subset(flag_subset$uniqueID,nchar(as.character(flag_subset$Sequence)) < seq_min)

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Min_Len"

  }

  #Here I get entries which don't have a flag
  flag_subset<-as.data.frame(subset(clean_data_table, clean_data_table$Flags == "-"), drop=FALSE)

  if(nrow(flag_subset)>0){

    #Here I flag entries which have sequences either above the desired size
    flag_subset<-subset(flag_subset$uniqueID,nchar(as.character(flag_subset$Sequence)) > seq_max )

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Max_Len"

  }

  #make all of the sequence characters caps
  clean_data_table$Sequence<-toupper(clean_data_table$Sequence)

  #remove all gaps from the sequences
  clean_data_table$Sequence<- gsub( "-", "", as.character(clean_data_table$Sequence))

  #Here I get entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  if(nrow(flag_subset>0)){

    #Flagging records with sequences with non IUPAC characters
    flag_subset<-subset(flag_subset$uniqueID, grepl("[^ACGTRYSWKMBDHVN-]", flag_subset$Sequence))

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Non-IUPAC"

  }

  #Here I get entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  if(nrow(flag_subset>0)){

    #Here I flag entries which don't have molecular marker information
    flag_subset<-subset(flag_subset$uniqueID, flag_subset$Gene=="")

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "No_Marker"

  }


  #Here I get entries which don't have a flag
  flag_subset<-subset(clean_data_table,clean_data_table$Flags == "-")

  if(nrow(flag_subset>0)){

    #Here I flag entries which don't have molecular marker information
    flag_subset<-subset(flag_subset$uniqueID, flag_subset$Genus!=target_genus)

    #Adding the results of the above check to the clean_data_table Flags column
    clean_data_table$Flags[clean_data_table$uniqueID %in% c(flag_subset)]<- "Non-Target"

  }

  #make all of the gene identifiers upper case
  clean_data_table$Gene<-toupper(clean_data_table$Gene)
  #remove all gaps from the gene names
  clean_data_table$Gene<- gsub( "[/]", "", as.character(clean_data_table$Gene))

  return(clean_data_table)
}
