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
    Main_Seq_file_data_frame<- Main_Seq_file_data_frame[!duplicated(Main_Seq_file_data_frame$Header),]
    
    #For each genus in the data file
    unique_genera_list <- unique(Main_Seq_file_data_frame$Genus)
    
    for(unique_genera in 1:length(unique_genera_list)){
      
      #Get the dataframe with the records from the loop genus
      Seq_file_data_frame <- Main_Seq_file_data_frame[Main_Seq_file_data_frame$Genus == unique_genera_list[unique_genera],]
      
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
          
          #This is adding the row names back to the matrix as the above line removed it.
          row.names(dist_matrix)<-colnames(dist_matrix)
          
          #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
          no_outliers_dist_matrix <- dist_matrix[!(rownames(dist_matrix) %in% seq_to_remove),]
          no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove)]
          
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
          
        }else{#end of the species more than 2 records check
          
          #Add a column with NA for the Species_Outlier spot
          log_df$Species_Outlier <- "NA"
        }
        #**************************** Setting up the Barcode Gap analysis section *********************************
        
        if(nrow(Seq_file_data_frame)>2){
          
          #Remove outliers from the calculated distance matrix now using all of the seq_to_remove obtained
          no_outliers_dist_matrix <- dist_matrix[ !(rownames(dist_matrix) %in% seq_to_remove),]
          no_outliers_dist_matrix <- no_outliers_dist_matrix[ , !(colnames(no_outliers_dist_matrix) %in% seq_to_remove)]
          
          #Remove the identified outliers from the main Seq file now using all of the seq_to_remove obtained
          barcode_gap_data_frame<-Seq_file_data_frame[!(Seq_file_data_frame$Header %in% seq_to_remove),]
          
          #Sort the data frame by Species
          barcode_gap_data_frame<- barcode_gap_data_frame[order(barcode_gap_data_frame$Species),]
          
          #Get a unique species list
          Species<-unique(barcode_gap_data_frame$Species)
          
          #Add barcode gap reporting columns to the species output reporting table log_df
          barcode_gap_columns<-c("Target_Records",
                                 "Target_Haplotypes",
                                 "Genus_Records",
                                 "Genus_Haplotypes",
                                 "Num_Overlap_Haploptype",
                                 "Overlap_Haploptype",
                                 "Intraspecific",
                                 "Interspecific",
                                 "Barcode_Gap",
                                 "Num_Barcode_Gap_Overlap_Records",
                                 "Barcode_Gap_Overlap_Records",
                                 "Num_Barcode_Gap_Overlap_Taxa",
                                 "Barcode_Gap_Overlap_Taxa", "Barcode_Gap_Result",
                                 "p_x",
                                 "q_x",
                                 "p_x_prime_NN",
                                 "q_x_prime_NN")
          log_df[,barcode_gap_columns] <- "-"
          
        }
        
        #*********************** Begin looping through the species in the genus ***************************
        
        if(nrow(Seq_file_data_frame)>2){
          
          for (species_list_counter in 1:length(Species)){
            
            #Get the records for the loop species
            loop_species_records<-barcode_gap_data_frame[barcode_gap_data_frame$Species==Species[species_list_counter],]
            
            #for species with more than one records so that the within is able to be calculated
            if (nrow(loop_species_records)>1 && length(Species)>1){
              
              
              
              
              
              
              ### Step [1] ###
              # Generate one vector that consists of all interspecific differences across all species pairs
              res <- outer(Species, Species, "==")
              # Since matrix is symmetric, replace upper triangle with NA and then subset to get lower triangular matrix
              no_outliers_dist_matrix[upper.tri(no_outliers_dist_matrix, diag = FALSE)] <- NA
              inter <- na.omit(no_outliers_dist_matrix[!res])
              ###############
              
              ### Step [2] ###
              #Get the rows of the target species from the dist matrix and then get the columns from the selected columns
              loop_species_dist_matrix <- no_outliers_dist_matrix[(rownames(no_outliers_dist_matrix) %in% loop_species_records$Header),]
              loop_species_dist_matrix_within <- loop_species_dist_matrix[,(colnames(loop_species_dist_matrix) %in% loop_species_records$Header)]
              intra <- na.omit(as.vector(loop_species_dist_matrix_within))
              ################
             
              
              #Now get comparisons between the loop species and all other records
              loop_species_dist_matrix_between <- loop_species_dist_matrix[,!(colnames(loop_species_dist_matrix) %in% loop_species_records$Header), drop=FALSE]
              
              
              ##### Steps [3] and [4] #####
              # split interspecific distance matrix by species
              # is.na(no_outliers_dist_matrix_lower) <- !res
              splt <- split(no_outliers_dist_matrix, sub("(?:(.*)\\|){2}(\\w+)\\|(\\w+)\\|.*?$", "\\1-\\2", colnames(no_outliers_dist_matrix)))
              
              # compute proportional overlap for nearest neighbours using mean distance
              splt1 <- lapply(splt, mean) 

              # convert to dataframe
              x <- as.data.frame(unlist(splt1))
              colnames(x) <- "Mean"

              # Pair each species with its nearest neighbour
              d <- data.frame(`diag<-`(as.matrix(dist(x$Mean)), Inf))
              ids <- unlist(Map(which.min, d)) # pair focal species index with its nearest neighbour index
              Neighbour <- x$Mean[ids]
              x <- data.frame(names(splt), x$Mean, Neighbour)
              names(x)[1] <- "Species"
              names(x)[2] <- "Mean Intraspecific Distance"
              x[, 3] <- x$Species[ids]

              splt2 <- splt[c(t(x[, c("Species", "Neighbour")]))] # rearrange list of distances so that focal species and nearest neighbours occur together
              ##########
              
              ##### Step [5] #####
         
              # compute proportional overlap of intraspecific and interspecific distributions
              
              # p_x is overlap of intra with inter
              p_x <- length(which(intra >= min(inter))) / length(intra)
              
              # q_x is overlap of inter with intra
              q_x <- length(which(inter <= max(intra))) / length(inter)
              
              # # loop through the splt2 list to compute p' and q'
              # # target species are in odd positions, nearest neighbours are in even positions
              # for (i in 1:length(splt2)) {
              #   for (j in 1:length(splt2)) {
              #     if ((i %% 2 == 1) && (j %% 2 == 0)) {
              #       p_x_prime_NN <- length(which(splt2[[i]] >= min(splt2[[j]]))) / length(splt2[[i]])
              #       q_x_prime_NN <- length(which(splt2[[j]] <= max(splt2[[i]]))) / length(splt2[[j]])
              #     }
              #   }
              # }
              
              ##########
              
              
              
              
              
              
              #Get the number of records for the target species
              loop_species_target<-nrow(loop_species_dist_matrix_within)
              
              #Get the unique target sequences for the haplotype reporting
              loop_target_haplotypes <- length(unique(barcode_gap_data_frame[barcode_gap_data_frame$Header %in% loop_species_records$Header, "Sequence"]))
              
              #There was clearly more than one record as we made it here past the more than one species record if so that means that the two records or more were
              #dereplicated. So, then that would mean that the between species distance will be 0 and the only record would have this as the comparison
              #Create the variable to report
              
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
              # p_x <- "-"
              # q_x <- "-"
              # p_x_prime_NN <- "-"
              # q_x_prime_NN <- "-"
              
            }#closing the if else checking if there is more than one species record
            
            #add the results of the species barcode gap within or intraspecific to the log_df
            log_df$Target_Records[log_df$Species %in% Species[species_list_counter] ]<-loop_species_target
            
            #Add the results of the target species haplotypes to the log df
            log_df$Target_Haplotypes[log_df$Species %in% Species[species_list_counter] ]<-loop_target_haplotypes
            
            #add the results of p_x
            log_df$p_x[log_df$Species %in% Species[species_list_counter] ]<- p_x
            
            #add the results of q_x
            log_df$q_x[log_df$Species %in% Species[species_list_counter] ]<- q_x
            
            #add the results of p_x_prime_NN
            #log_df$p_x_prime_NN[log_df$Species %in% Species[species_list_counter] ]<- p_x_prime_NN
            
            #add the results of q_x_prime_NN
            #log_df$q_x_prime_NN[log_df$Species %in% Species[species_list_counter] ]<- q_x_prime_NN
            
            
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
            
           
            ### DNA barcode gap overlap plots ###
            
            # convert to log-10 scale
            # since p and q can be 0, plotting on log10 scale allows easier visualization
            p_x <- as.numeric(log_df$p_x)
            q_x <- as.numeric(log_df$q_x)
            # p_x_prime_NN <- as.numeric(log_df$p_x_prime_NN)
            # q_x_prime_NN <- as.numeric(log_df$q_x_prime_NN)
            
            # replace infinite values with -5
            df_pq <- data.frame(log10(p_x), log10(q_x))
            df_pq <- apply(df_pq, 2, function(x) replace(x, is.infinite(x), -5)) # replace Inf with finite value
            df_pq <- as.data.frame(df_pq)
            
            # df_pq_prime_NN <- data.frame(log10(p_x_prime_NN), log10(q_x_prime_NN))
            # df_pq_prime_NN <- apply(df_pq_prime_NN, 2, function(x) replace(x, is.infinite(x), -5))
            # df_pq_prime_NN <- as.data.frame(df_pq_prime_NN)
            
            # Plot the results
            p <- ggplot(df_pq, aes(x = log10.p_x., y = log10.q_x.)) + geom_point(colour = "blue") +
              labs(x = expression(log[10](p)), y = expression(log[10](q)))
            
            #save plot to file without using ggsave
            png(paste0(Work_loc,"/",file_name[h],"_pq.png"))
            print(p)
            dev.off()
            
            # p <- ggplot(df_pq_prime_NN, aes(x = log10.p_x_prime_NN., y = log10.q_x_prime_NN.)) + geom_point(colour = "blue") +
            #   labs(x = expression(log[10](p*"'")), y = expression(log[10](q*"'")))
            
            # png(paste0(Work_loc,"/",file_name[h],"_pq_prime_NN.png"))
            # print(p)
            # dev.off()
            
          }
          
        } #End of species in genus section
        
      }#End of the if more than two records in the genus loop dataset
      
    }#End of the genus loop
    
  } #end of file loop
  
  print(paste("Start time... ",start_time," and end time... ",Sys.time()))

