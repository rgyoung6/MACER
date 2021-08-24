# Molecular Acquisition, Cleaning, and Evaluation in R (MACER) information and example.

# Description:

This repository contains the MACER package located at rgyoung6/MACER. The Molecular Acquisition, Cleaning, and Evaluation in R (MACER) is a tool to assemble, align, trim, and evaluate molecular sequence datasets from BOLD and GenBank.

# Functions:
**auto_seq_download()**<br/>
**create_fastas()**<br/>
**align_to_ref()**<br/>
**barcode_clean()**<br/>

# Installation 

MACER can be installed three ways.

## 1. Install from CRAN

install.packages('MACER’)

## 2. Install via GitHub
Run the following commands in your R terminal...<br/>

> install.packages("devtools")<br/>
> library(devtools)<br/>
> devtools::install_github("rgyoung6/MACER")<br/>

**Note:** the first command to install the "devtools" may not be necessary if already installed.<br/>

## 3. Install through download from GitHub.
Navigate to the [MACER](https://github.com/rgyoung6/MACER) GitHub page. Download the files associated with this page to your local computer and place them somewhere in the main file folder named MACER. Then run the following command pointing to that location on your local computer by replacing the HERE with the path in the below command...<br/>
> library("MACER", lib.loc="HERE")<br/>

# Function Descriptions:

## auto_seq_download()
This function takes a list of genera, as supplied by the user, and searches and downloads molecular sequence data from BOLD and Genbank.

## create_fastas()
This function takes an output file from the auto_seq_download() function and a user supplied table of taxa and marker names desired and outputs fasta files with target records.

## align_to_ref()
This function takes a FASTA file with target sequences and aligns them against a reference sequence submitted to the program. The output is an aligned fasta file that is trimmed to the length of the reference sequence. Sequences without full coverage (records having sequences with leading or trailing gaps) are removed. Finally, internal gaps are removed from the sequence based on the submitted multiple sequence alignment (MSA) percent coverage of the character position as provided in the pigl argument supplied by the user.

## barcode_clean()
This function takes an input fasta file and removes genus level outliers and species outliers based on the 1.5 x greater than the interquartile range both at the genus and species level. It also, if selected, checks the sequence using amino acid translation and has the option to eliminate sequences that have non-IUPAC codes. Finally, the program calculates the barcode gap for the species in the submitted dataset.

# Running MACER 

## NOTE: When running MACER scripts the paths for the files selected cannot have whitespace! Any file with whitespace in the file naming may cause the program to not run and end in an error, or may result in unexpected output.

# Package pipeline example:
An example walk through of the MACER functions is available in the associated GitHub repository rgyoung6/MACER_example. 

# Package Function Details

## Download – auto_seq_download()

### Input 
File with list of genera you want to download. Place these genera in a file in a single column.

### Arguments
BOLD_database – TRUE is to include, FALSE is to exclude; default TRUE
NCBI_database – TRUE is to include, FALSE is to exclude; default TRUE
search_str – NULL uses the default string, anything other than NULL then that string will be used for the GenBank search; default NULL
Default string…
(genus[ORGN]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])
Note: When using a custom search string for NCBI only a single genus at a time can be used.
input_file –  FALSE prompts the user to indicate the location of the input file through point and click prompts, anything other than FALSE then the string supplied will be used for the location; default FALSE
output_file –  FALSE prompts the user to indicate the location of the output file through point and click prompts, anything other than FALSE then the string supplied will be used for the location; default FALSE

### Output
One main folder containing three other folders.
Main folder - Seq_auto_dl_TTTTTT_MMM_DD
Three subfolders
BOLD - Contains a file for every genus downloaded with the raw data from the BOLD system.
NCBI - Contains a file for every genus downloaded with the raw data from GenBank.
Total_tables - Contains files for the running of the script these include.
A_Summary.txt - This file contains information about the downloads.
A_Total_Table.tsv – This is a file with a single table (tab delimited) containing the accumulated data obtained for all genera searched.

Note: the A_Total_Table.tsv file contains all of the records obtained and the final column is a flagged file of the results of the auto_seq_download() function

### Dependencies
rentrez used to access and download sequences from NCBI’s GenBank

## Table to Fasta – create_fastas()

### Input 
File with list of genera with the molecular markers names below the taxa. The information to create this parameters file can be obtained from A_Summary.txt file from the download script results. For example see below…
Tamias<t/> Tamias
CYTB<t/> COI-5P
CYTOCHROMEB<t/>	CYTOCHROMECOXIDASESUBUNIT1
CYTOCHROME-B<t/>	CYTOCHROMECOXIDASESUBUNITI

### Arguments
data_file – FALSE prompts the user to indicate the location of the data file in the format of the auto_seq_download output, anything other than FALSE then the string supplied will be used for the location; default FASLE
input_file – FALSE prompts the user to indicate the location of the input file used to select through point and click prompts, anything other than FALSE then the string supplied will be used for the location; default FALSE
output_folder – FALSE prompts the user to indicate the location of the output file through point and click prompts, anything other than FALSE then the string supplied will be used for the location; default FALSE

no_marker – If set to TRUE then will include records filtered out due to no marker data. Default is FALSE to not include records with no marker data.
no_taxa – If set to TRUE then will include records filtered out due to no taxa data. Default is FALSE to not include records with no taxa data.
no_seq – If set to TRUE then will include records filtered out due to no sequence data. Default is FALSE to not include records with no sequence data.
name_issue – If set to TRUE then will include records filtered out due to genus and species names with more than two terms. Default is FALSE to not include records with taxonomic naming issues.
taxa_digits – If set to TRUE then will include records filtered out due to genus or species names containing digits. Default is FALSE to not include records with digits in the taxonomic naming.
taxa_punct – If set to TRUE then will include records filtered out due to the presence of punctuation in the genus or species names. Default is FALSE to not include records with punctuation in the taxonomic naming.
wrong_taxa – If set to TRUE then will include records filtered out due to the incorrect genera based on the list of taxa initially submitted to the download program. Default is FALSE to not include records of non-target taxa.

### Output
This script outputs a fasta file of sequences for each column in the submitted parameters file. These files are named with the genera of interest and the first marker name in the column of the parameters file. These files are located in the folder where the Total_tables.tsv file is located.

### Dependencies
N/A

## Align - align_to_ref()

### Input 
A file folder location with the fasta files that need to be aligned and trimmed using the supplied reference sequence. Please note that any and all fasta files (named *.fas) in this folder will be analyzed.
A reference sequence file with a sequence or MSA with all sequences having the same length. 
The location of the MAFFT executable file (https://mafft.cbrc.jp/alignment/software/)

### Arguments
data_folder – This variable can be used to provide a location for the file containing all of the fasta files wanting to be aligned. The default value is set to FALSE where the program will prompt the user to select the folder through point-and-click.
ref_seq_file – This variable can be used to provide a location for the reference sequence file. The default value is set to FALSE where the program will prompt the user to select the folder through point-and-click.
MAFFT_loc – This variable can be used to provide a location for the MAFFT program. The default value is set to FALSE where the program will prompt the user to select the folder through point-and-click.
output_file – This variable can be used to set the location of the output files from the program. The default value is set to FALSE where the program will place the output files in the same location as the target files.

pigl – This is the percent internal gap loop argument. This provides a percent that will remove records causing internal gaps is more than the percent value assigned to this argument is reached. If this value is set to 0 then internal gaps are not removed. The default for this value is 0.95.
op – This is the gap opening penalty for the use of MAFFT. The higher the value the larger penalty in the alignment. The default for this value is set to 10. The default value in the MAFFT program is 1.53. For alignment of highly conserved regions where no gaps are expected this should be set to a much higher number and 10 is recommended for barcode regions like the COI-5P.

### Output
In the submitted file folder location there will be a log file titled “MAFFT_log”.
The sequence output files from this script are placed into two subfolders. These folders are in the submitted file location where the fasta files of interest are located. The two folders created are MAFFT and MAFFT_trimmed. In the MAFFT folder there will be files with name of the files in the submitted file folder appended with “_MAFFT”. 
The MAFFT_trimmed file will contain files with the same naming convention as the files in the submitted folder and appended with “MAFFT_trimmed”.

### Dependencies
While not an R dependency, the use of the MAFFT program as installed on the local computer is required. 

## Clean - barcode_clean()

### Input 
A file folder with one or more fasta files of interest

### Arguments
AA_code – This is the amino acid translation matrix used to check the sequences for stop codons. The following codes are available. The Invertebrate matrix, 2, is the default.
  - std (1 in ape) is standard code
  - vert (2 in ape) is vertebrate mitochondrial
  - invert (5 in ape) is invertebrate mitochondrial
  - F skips the AA clean section
AGCT_only – This only keeps sequences with AGCT exclusively, not IUPAC characters.
  - TRUE is on 
  - FALSE is accepting all IUPAC characters

### Output
A single log file for the running of the function with the name A_Clean_File_YYYY-DD-TTTTTTTT. The function will also output three files for each fasta file submitted. The first is the distance matrix that was calculated and used to assess the DNA barcode gaps. This file is named the same as the input file with ‘_dist_table.dat” appended to the end of the name. The second file is the total data table file which provides a table of all submitted records for each data set accompanied with the results from each section of the analysis. This file is named the same as the input fasta with “_data_table.dat” appended to the end, Finally, a fasta file with all outliers and flagged records removed is generated for each input fasta file. This output file is named the same as the input fasta with “_no_outlier.fas” appended to the end. The flags that are possible are non_AGCT, Stop_Codon, Genus_Outlier, Species_Outlier, and '-'.

### Dependencies
ape is required for distance matrix construction.
