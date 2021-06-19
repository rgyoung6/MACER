
# Molecular Acquisition, Cleaning, and Evaluation in R (MACER) information and example.

# Description:

This repository contains the example files and outputs from the example for the MACER package located at rgyoung6/MACER. The Molecular Acquisition, Cleaning, and Evaluation in R (MACER) is a tool to assemble, align, trim, and evaluate molecular sequence datasets from BOLD and GenBank.

# Functions:
**auto_seq_download()**<br/>
**create_fastas()**<br/>
**align_to_ref()**<br/>
**barcode_clean()**<br/>

# Installation 

MACER is not currently on Comprehensive R Archive Network (CRAN) and so you have two ways to use the package functions.

## 1. Install via GitHub
Run the following commands in your R terminal...<br/>

> install.packages("devtools")<br/>
> library(devtools)<br/>
> devtools::install_github("rgyoung6/MACER")<br/>

**Note:** the first command to install the "devtools" may not be necessary if already installed.<br/>

## 2. Install through download from GitHub.
Navigate to the [MACER](https://github.com/rgyoung6/MACER) GitHub page. Download the files associated with this page to your local computer and place them somewhere in the main file folder named MACER. Then run the following command pointing to that location on your local computer by replacing the HERE with the path in the below command...<br/>
> library("MACER", lib.loc="HERE")<br/>

# Function Descriptions:

## auto_seq_download()
This function takes a list of genera, as supplied by the user, and searches and downloads molecular sequence data from BOLD and Genbank.

## create_fastas()
This function takes an input fasta file and removes genus level outliers and species outliers based on the 1.5 x greater than the interquartile range. It also, if selected, checks the sequence using amino acid translation and has the option to eliminate sequences that have non-IUPAC codes. Finally, the program calculates the barcode gap for the species in the submitted dataset.

## align_to_ref()
This function takes a FASTA file with target sequences and aligns them against a reference sequence submitted to the program. The output is an aligned fasta file that is trimmed to the length of the reference sequence. Sequences without full coverage (records having sequences with leading or trailing gaps) are removed. Records with characters other than IUPAC are also removed. Finally, internal gaps are removed from the sequence based on the submitted multiple sequence alignment (MSA) percent coverage of the character position as provided in the pigl argument supplied by the user.

## barcode_clean()
This function takes an input fasta file and removes genus level outliers and species outliers based on the 1.5 x greater than the interquartile range. It also, if selected, checks the sequence using amino acid translation and has the option to eliminate sequences that have non-IUPAC codes. Finally, the program calculates the barcode gap for the species in the submitted dataset.

# Running MACER 

## NOTE: When running MACER scripts the paths for the files selected cannot have whitespace! Any file with whitespace in the file naming may cause the program to not run and end in an error, or may result in unexpected output.

# Package pipeline example:
<<<<<<< HEAD
An example walk through of the MACER functions is available in the associated GitHub repository rgyoung6/MACER_example. 

# Package Function Details

## Download – auto_seq_download()

### Input 
File with list of genera you want to download. Place these genera in a file in a single column.

### Arguments
BOLD_database – 1 is to include, 0 is to exclude; default 1
NCBI_database – 1 is to include, 0 is to exclude; default 1
search_str – N uses the default string, anything other than N then that string will be used for the GenBank search; default N.
Default string…
(genus[ORGN] OR genus[ALL]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])
Note: When using a custom search string for NCBI only a single genus at a time can be used.

### Output
One main folder containing three other folders.
Main folder - Seq_auto_dl_######_MMM_DD
Three subfolders
BOLD - Contains a file for every genus downloaded with the raw data from the BOLD system.
NCBI - Contains a file for every genus downloaded with the raw data from GenBank.
Total_tables - Contains files for the running of the script these include.
A_Summary.txt - This file contains information about the downloads.
A_Total_Table.txt – This is a file with a single table (tab delimited) containing the accumulated data obtained for all genera searched.

Note: the A_Total_Table.txt file contains all of the records obtained and the final column is a flagged file of the sults of the auto_seq_download() function

### Dependencies
rentrez ver. 1.2.2 used to access and download sequences from NCBI’s GenBank

## Table to Fasta – create_fastas()

### Input 
File with list of genera with the molecular markers names below the taxa. The information to create this parameters file can be obtained from A_Summary.txt file from the download script results. For example see below…
Tamias	Tamias
CYTB	COI-5P
CYTOCHROMEB	CYTOCHROMECOXIDASESUBUNIT1
CYTOCHROME-B	CYTOCHROMECOXIDASESUBUNITI

### Arguments
no_marker – If set to 1 then will include records filtered out due to no marker data. Default is 0 to not include records with no marker data.
no_taxa – If set to 1 then will include records filtered out due to no taxa data. Default is 0 to not include records with no taxa data.
no_seq – If set to 1 then will include records filtered out due to no sequence data. Default is 0 to not include records with no sequence data.
name_issue – If set to 1 then will include records filtered out due to genus and species names with more than two terms. Default is 0 to not include records with taxonomic naming issues.
taxa_digits – If set to 1 then will include records filtered out due to genus or species names containing digits. Default is 0 to not include records with digits in the taxonomic naming.
taxa_punct – If set to 1 then will include records filtered out due to the presence of punctuation in the genus or species names. Default is 0 to not include records with punctuation in the taxonomic naming.
wrong_taxa – If set to 1 then will include records filtered out due to the incorrect genera based on the list of taxa initially submitted to the download program. Default is 0 to not include records of non-target taxa.

### Output
This script outputs a fasta file of sequences for each column in the submitted parameters file. These files are named with the genera of interest and the first marker name in the column of the parameters file. These files are located in the folder where the Total_tables.dat file is located.

### Dependencies
N/A

## Align - align_to_ref()

### Input 
A file folder location with the fasta files that need to be aligned and trimmed using the supplied reference sequence. Please note that any and all fasta files (named *.fas) in this folder will be analyzed.
A reference sequence file with a sequence or MSA with all sequences having the same length. 
The location of the MAFFT executable file (https://mafft.cbrc.jp/alignment/software/)

### Arguments
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
  - 1 is standard code
  - 2 is vertebrate mitochondrial
  - 5 is invert mitochondrial
  - 0 skips the AA clean section
AGCT_only – This only keeps sequences with AGCT exclusively, not IUPAC characters.
  - 1 is on 
  - 0 is accepting all IUPAC characters

### Output
A single log file for the running of the function with the name A_Clean_File_########## where the numbers represent the date and time of the run. The function will also output three files for each fasta file submitted. The first is the distance matrix that was calculated and used to assess the DNA barcode gaps. This file is named the same as the input file with ‘_dist_table.dat” appended to the end of the name. The second file is the total data table file which provides a table of all submitted records for each data set accompanied with the results from each section of the analysis. This file is named the same as the input fasta with “_data_table.dat” appended to the end, Finally, a fasta file with all outliers and flagged records removed is generated for each input fasta file. This output file is named the same as the input fasta with “_no_outlier.fas” appended to the end. The flags that are possible are non_AGCT, Stop_Codon, Genus_Outlier, Species_Outlier, and '-'.

### Dependencies
ape 5.4-1 is required for distance matrix construction.
=======
The following example will walk through the use of the four main function in the MACER package. To illustrate the use of the package this example will use chipmunks as the example. The three genera, Eutamias, Neotamias, Tamias were used to search both BOLD and GenBank. Although Eutamias and Neotamias are not currently accepted taxonomic genera, these were retained to evaluate potential data with older taxonomic naming conventions. All associated files necessary to run this example, and the outputs from the running of this example are included in the GitHub MACER repository and the 'Example' file folder.

## Download target records using auto_seq_download()
Using the file Chipmunk.txt as the taxa for this example. This file contains a list of genera for Chipmunks, Eutamias, Neotamias, and Tamias. Run the auto_seq_download() function and when prompted hit enter and then point to the Chipmunk.txt file. The program will run and populate files in the same directory with the download results. While there are three potential arguments for this function this example will not use any. The reason arguments were not included is because the defaults of including BOLD and GenBank searches and using the default search algorithm were desired. The following is the output on the R screen for this example.

	> auto_seq_download()
	Choose the file with the genera of interest to download. Hit enter key to continue...
	"C:\\A_MACER\\Chipmunk.dat"
	"Starting time - 2021-06-16 14:51:20 - Eutamias"
	"Attempting BOLD Download - 2021-06-16 14:51:20 - Eutamias"
	"BOLD Download Error: failed to populate data table or there is no data"
	"Attempting NCBI Download - Eutamias"
	"Here is the value of the search string...(Eutamias[ORGN] OR Eutamias[ALL]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])"
	"Downloading 1000 Eutamias NCBI records starting at...1 of 3248 - 2021-06-16 14:51:22"
	"Downloading 1000 Eutamias NCBI records starting at...1001 of 3248 - 2021-06-16 14:52:05"
	"Downloading 1000 Eutamias NCBI records starting at...2001 of 3248 - 2021-06-16 14:53:19"
	"Downloading 1000 Eutamias NCBI records starting at...3001 of 3248 - 2021-06-16 14:54:28"
	"Download Complete: cleaning NCBI data table - Eutamias"
	"Starting time - 2021-06-16 14:54:39 - Tamias"
	"Attempting BOLD Download - 2021-06-16 14:54:39 - Tamias"
	"Download Complete: cleaning BOLD data table - Tamias"
	"Attempting NCBI Download - Tamias"
	"Here is the value of the search string...(Tamias[ORGN] OR Tamias[ALL]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])"
	"Downloading 1000 Tamias NCBI records starting at...1 of 3459 - 2021-06-16 14:54:46"
	"Downloading 1000 Tamias NCBI records starting at...1001 of 3459 - 2021-06-16 14:55:20"
	"Downloading 1000 Tamias NCBI records starting at...2001 of 3459 - 2021-06-16 14:56:31"
	"Downloading 1000 Tamias NCBI records starting at...3001 of 3459 - 2021-06-16 14:57:58"
	"Download Complete: cleaning NCBI data table - Tamias"
	"Starting time - 2021-06-16 14:58:21 - Neotamias"
	"Attempting BOLD Download - 2021-06-16 14:58:21 - Neotamias"
	"BOLD Download Error: failed to populate data table or there is no data"
	"Attempting NCBI Download - Neotamias"
	"Here is the value of the search string...(Neotamias[ORGN] OR Neotamias[ALL]) NOT (shotgun[ALL] OR genome[ALL] OR assembled[ALL] OR microsatellite[ALL])"
	"Downloading 1000 Neotamias NCBI records starting at...1 of 3245 - 2021-06-16 14:58:23"
	"Downloading 1000 Neotamias NCBI records starting at...1001 of 3245 - 2021-06-16 14:59:01"
	"Downloading 1000 Neotamias NCBI records starting at...2001 of 3245 - 2021-06-16 15:00:09"
	"Downloading 1000 Neotamias NCBI records starting at...3001 of 3245 - 2021-06-16 15:01:31"
	"Download Complete: cleaning NCBI data table - Neotamias"

The results obtained from the download script indicates that there were records present in the NCBI-GenBank database with each of the genera names. However, there were no records in the BOLD database for Eutamias and Neotamias. The output file A_Summary.txt from the auto_seq_download() function will contain information on the downloaded records. This file will be located in the Seq_auto_dl_######_MMM_DD file folder and the subfolder Total_Tables. For this examples case the file folder was named Seq_auto_dl_063034_Jun_07. The following is the contents of this file from this example.

	Starting time - 2021-06-07 06:30:34 - Eutamias
	Attempting BOLD Download - 2021-06-07 06:30:34 - Eutamias
	Attempting NCBI Download - 2021-06-07 06:30:35 - Eutamias
	Attempting NCBI Download - 2021-06-07 06:30:35 - Eutamias
	Number of records - 3355
	Species - 	amoenus,minimus,sibiricus,striatus,palmeri,monacensis,lusitaniae,senex,speciosus,alpinus,ruficaudus,townsendii,umbrinus,sonomae,siskiyou,rufus,quadrivittatus,quadrimaculatus,panamintinus,ochrogenys,obscurus,merriami,durangae,dorsalis,cinereicollis,canipes,bulleri,muris,-
	Molecular markers - CYTOCHROMEB,CYTOCHROME-B,CYTOCHROMECOXIDASESUBUNITII,16SRIBOSOMALRNA,,SIRTUIN6,CYTOCHROMEOXIDASESUBUNITI,VONWILLEBRANDFACTOR,RECOMBINATIONACTIVATINGPROTEIN1,GROWTHHORMONERECEPTOR,ENAMELIN,APOLIPOPROTEINB,ALPHA2BADRENERGICRECEPTOR,INTERPHOTORECEPTORRETINOIDBINDINGPROTEIN,EDG1,12SRIBOSOMALRNA,TYROSINASE,RECOMBINATIONACTIVATINGPROTEIN2,PREPRONOCICEPTIN,PHOSPHOLIPASECBETA4,CAMPRESPONSIVEELEMENTMODERATOR,CANNABINOIDRECEPTOR1,BMI1,BRAIN-DERIVEDNEUROTROPHICFACTOR,ATP7A,AMYLOIDBETAPRECURSORPROTEIN,ADENOSINEA3RECEPTOR,BETA-2ADRENERGICRECEPTOR,CITRATESYNTHASE,CYTOCHROMEOXIDASESUBUNIT1,NUCLEARRECEPTORSUBFAMILY0GROUPBMEMBER2,ZINCFINGERPROTEINZFX,C-MYC,C-MYCPROTEIN,ACROSIN,ACIDPHOSPHATASE5,TRNA-PHE,PEROXISOMEPROLIFERATOR-ACTIVATEDRECEPTORGAMMA,RAG1PROTEIN,THYROGLOBULIN,BETASPECTRIN,PRKC1,THYROTROPINBETASUBUNIT,MASTICATORYSUPERFASTMYOSINHEAVYCHAIN,MASTICATORYSUPERFASTMYOSINLIGHTCHAIN2,EMBRYONICATRIALMYOSINLIGHTCHAIN1,INTERPHOTORECEPTORBINDINGPROTEIN,ENVELOPEPROTEINSYNCYTIN-MAR1,CYTOCHROMEBOXIDASE,CYTOCHROMECOXIDASESUBUNIT2,LACTATEDEHYDROGENASEC,HEMOGLOBINBETA,HEMOGLOBINALPHA,SMCYPROTEIN,BREASTCANCERSUSCEPTIBILITY1,BREASTANDOVARIANCANCERSUSCEPTIBILITY1,DENTINMATRIXPROTEIN1,CYTOCHROMECOXIDASESUBUNITI,ZONAPELLUCIDAGLYCOPROTEIN2,ZONADHESIN,TTN,GPROTEINBETASUBUNIT5SHORTVARIANT,CONEPHOSPHODIESTERASEBETASUBUNIT,ORFII,REGULATOROFGPROTEINSIGNALINGRGS9-1,HIBERNATIONSPECIFICPROTEIN27,FBN1,BCHE,HP-55,ALPHA1-ANTITRYPSIN-LIKEPROTEIN,INTERPHOTORECEPTERRETINOIDBINDINGPROTEIN,RAG1,HIBERNATIONSPECIFICPROTEIN-25,MGF,SMALLSUBUNITRIBOSOMALRNA,HP-20,HP-25,HEPATOCYTENUCLEARFACTOR4,NADHDEHYDROGENASESUBUNIT1,HP-27,TRANSCRIPTIONFACTORSP1,GLYCERALDEHYDE-3-PHOSPHATEDEHYDROGENASE,HEATSHOCKFACTOR1,HEATSHOCK70KDAPROTEIN1A,SERUMALBUMINPREPROTEIN,SERUMALBUMINPREPROPROTEIN
	Ending time - 2021-06-07 06:35:17 - Eutamias
	Starting time - 2021-06-07 06:35:17 - Tamias
	Attempting BOLD Download - 2021-06-07 06:35:17 - Tamias
	Download Complete: cleaning BOLD data table - 2021-06-07 06:35:30 - Tamias
	Attempting NCBI Download - 2021-06-07 06:35:30 - Tamias
	Attempting NCBI Download - 2021-06-07 06:35:30 - Tamias
	Number of records - 3741
	Species - sibiricus,striatus,amoenus,minimus,umbrinus,quadrivittatus,dorsalis,rufus,cinereicollis,canipes,ruficaudus,phagocytophilum,folkertsi,burgdorferi,lanei,bissettiae,niloticus,sapiens,acanthias,monax,palmeri,senex,speciosus,alpinus,townsendii,sonomae,siskiyou,quadrimaculatus,panamintinus,ochrogenys,obscurus,merriami,durangae,bulleri,erratica,parvum,hudsonicus,microti,tamias,ezoensis,hermsii,sciuricola,washoeensis,bovis,unidentified,-
	Molecular markers - COI-5P,ND3,ND4L,ND6,COII,ND2,ND5-0,COXIII,CYTB,ND1,ND4,CYTOCHROMEB,CYTOCHROME-B,CYTOCHROMECOXIDASESUBUNITII,HEATSHOCKPROTEIN,16SRIBOSOMALRNA,FLAGELLIN,,SIRTUIN6,CYTOCHROMEOXIDASESUBUNITI,VONWILLEBRANDFACTOR,RECOMBINATIONACTIVATINGPROTEIN1,GROWTHHORMONERECEPTOR,ENAMELIN,APOLIPOPROTEINB,ALPHA2BADRENERGICRECEPTOR,INTERPHOTORECEPTORRETINOIDBINDINGPROTEIN,EDG1,12SRIBOSOMALRNA,TYROSINASE,RECOMBINATIONACTIVATINGPROTEIN2,PREPRONOCICEPTIN,PHOSPHOLIPASECBETA4,CAMPRESPONSIVEELEMENTMODERATOR,CANNABINOIDRECEPTOR1,BMI1,BRAIN-DERIVEDNEUROTROPHICFACTOR,ATP7A,AMYLOIDBETAPRECURSORPROTEIN,ADENOSINEA3RECEPTOR,BETA-2ADRENERGICRECEPTOR,CYTOCHROMEOXIDASESUBUNIT1,NUCLEARRECEPTORSUBFAMILY0GROUPBMEMBER2,ZINCFINGERPROTEINZFX,C-MYC,C-MYCPROTEIN,ACROSIN,ACIDPHOSPHATASE5,TRNA-PHE,NADHDEHYDROGENASESUBUNIT1,PEROXISOMEPROLIFERATOR-ACTIVATEDRECEPTORGAMMA,ALPHA1-ANTITRYPSIN-LIKEPROTEIN,RAG1PROTEIN,THYROGLOBULIN,BETASPECTRIN,PRKC1,THYROTROPINBETASUBUNIT,ELONGATIONFACTOR1-ALPHA,CYTOCHROMECOXIDASESUBUNITI,MASTICATORYSUPERFASTMYOSINHEAVYCHAIN,MASTICATORYSUPERFASTMYOSINLIGHTCHAIN2,EMBRYONICATRIALMYOSINLIGHTCHAIN1,INTERPHOTORECEPTORBINDINGPROTEIN,CITRATESYNTHASE,POLYTHREONINE-RICHGLYCOPROTEIN,DIHYDROFOLATEREDUCTASE,18SRIBOSOMALRNA,HISTONEH3,CYTOCHROMECOXIDASESUBUNIT1,28SRIBOSOMALRNA,FACTORHBINDINGPROTEIN,VARIABLETICKPROTEIN,GLYCEROPHOSPHODIESTERPHOSPHODIESTERASE,DNAGYRASESUBUNITB,16S-23SRIBOSOMALRNAINTERGENICSPACER,ENVELOPEPROTEINSYNCYTIN-MAR1,CYTOCHROMEBOXIDASE,CYTOCHROMECOXIDASESUBUNIT2,LACTATEDEHYDROGENASEC,HEMOGLOBINBETA,HEMOGLOBINALPHA,HEATSHOCKPROTEINCOGNATE5,SIGNALRECOGNITIONPARTICLEPROTEIN54K,SMCYPROTEIN,BREASTCANCERSUSCEPTIBILITY1,BREASTANDOVARIANCANCERSUSCEPTIBILITY1,DENTINMATRIXPROTEIN1,ZONAPELLUCIDAGLYCOPROTEIN2,ZONADHESIN,TTN,OUTERSURFACEPROTEIN,RNAPOLYMERASEBETASUBUNIT,RIBOFLAVINSYNTHASE,60KDAHEATSHOCKPROTEIN,CELLDIVISIONPROTEIN,18SSMALLSUBUNITRIBOSOMALRNA,GPROTEINBETASUBUNIT5SHORTVARIANT,CONEPHOSPHODIESTERASEBETASUBUNIT,ORFII,REGULATOROFGPROTEINSIGNALINGRGS9-1,HIBERNATIONSPECIFICPROTEIN27,FBN1,BCHE,HP-55,INTERPHOTORECEPTERRETINOIDBINDINGPROTEIN,RAG1,HIBERNATIONSPECIFICPROTEIN-25,MGF,HP-20,SMALLSUBUNITRIBOSOMALRNA,HP-25,HEPATOCYTENUCLEARFACTOR4,HP-27,POLYMERASE,TRANSCRIPTIONFACTORSP1,GLYCERALDEHYDE-3-PHOSPHATEDEHYDROGENASE,HEATSHOCKFACTOR1,HEATSHOCK70KDAPROTEIN1A,SERUMALBUMINPREPROTEIN,SERUMALBUMINPREPROPROTEIN,RPOB,OUTERMEMBRANEPROTEINA,17-KDAGENUSSPECIFICANTIGEN,ACTIN,ENVELOPEPROTEIN,FT'55MSPROTEASEINHIBITOR'FT,FT'HP-55PROTEASEINHIBITOR'FT,FT'55RSPROTEASEINHIBITOR'FT,FT'55MMPROTEASEINHIBITOR'FT
	Ending time - 2021-06-07 06:40:01 - Tamias
	Starting time - 2021-06-07 06:40:01 - Neotamias
	Attempting BOLD Download - 2021-06-07 06:40:01 - Neotamias
	Attempting NCBI Download - 2021-06-07 06:40:01 - Neotamias
	Attempting NCBI Download - 2021-06-07 06:40:01 - Neotamias
	Number of records - 3352
	Species - amoenus,minimus,sibiricus,striatus,palmeri,senex,speciosus,alpinus,ruficaudus,townsendii,umbrinus,sonomae,siskiyou,rufus,quadrivittatus,quadrimaculatus,panamintinus,ochrogenys,obscurus,merriami,durangae,dorsalis,cinereicollis,canipes,bulleri,lateralis,-
	Molecular markers - CYTOCHROMEB,CYTOCHROME-B,CYTOCHROMECOXIDASESUBUNITII,16SRIBOSOMALRNA,,SIRTUIN6,CYTOCHROMEOXIDASESUBUNITI,VONWILLEBRANDFACTOR,RECOMBINATIONACTIVATINGPROTEIN1,GROWTHHORMONERECEPTOR,ENAMELIN,APOLIPOPROTEINB,ALPHA2BADRENERGICRECEPTOR,INTERPHOTORECEPTORRETINOIDBINDINGPROTEIN,EDG1,12SRIBOSOMALRNA,TYROSINASE,RECOMBINATIONACTIVATINGPROTEIN2,PREPRONOCICEPTIN,PHOSPHOLIPASECBETA4,CAMPRESPONSIVEELEMENTMODERATOR,CANNABINOIDRECEPTOR1,BMI1,BRAIN-DERIVEDNEUROTROPHICFACTOR,ATP7A,AMYLOIDBETAPRECURSORPROTEIN,ADENOSINEA3RECEPTOR,BETA-2ADRENERGICRECEPTOR,CYTOCHROMEOXIDASESUBUNIT1,NUCLEARRECEPTORSUBFAMILY0GROUPBMEMBER2,ZINCFINGERPROTEINZFX,C-MYC,C-MYCPROTEIN,ACROSIN,ACIDPHOSPHATASE5,TRNA-PHE,PEROXISOMEPROLIFERATOR-ACTIVATEDRECEPTORGAMMA,RAG1PROTEIN,THYROGLOBULIN,BETASPECTRIN,PRKC1,THYROTROPINBETASUBUNIT,MASTICATORYSUPERFASTMYOSINHEAVYCHAIN,MASTICATORYSUPERFASTMYOSINLIGHTCHAIN2,EMBRYONICATRIALMYOSINLIGHTCHAIN1,INTERPHOTORECEPTORBINDINGPROTEIN,ENVELOPEPROTEINSYNCYTIN-MAR1,CYTOCHROMEBOXIDASE,CYTOCHROMECOXIDASESUBUNIT2,LACTATEDEHYDROGENASEC,HEMOGLOBINBETA,HEMOGLOBINALPHA,SMCYPROTEIN,BREASTCANCERSUSCEPTIBILITY1,BREASTANDOVARIANCANCERSUSCEPTIBILITY1,DENTINMATRIXPROTEIN1,CYTOCHROMECOXIDASESUBUNITI,ZONAPELLUCIDAGLYCOPROTEIN2,ZONADHESIN,TTN,GPROTEINBETASUBUNIT5SHORTVARIANT,CONEPHOSPHODIESTERASEBETASUBUNIT,ORFII,REGULATOROFGPROTEINSIGNALINGRGS9-1,HIBERNATIONSPECIFICPROTEIN27,FBN1,BCHE,HP-55,ALPHA1-ANTITRYPSIN-LIKEPROTEIN,INTERPHOTORECEPTERRETINOIDBINDINGPROTEIN,RAG1,HIBERNATIONSPECIFICPROTEIN-25,MGF,HP-20,HP-25,HEPATOCYTENUCLEARFACTOR4,NADHDEHYDROGENASESUBUNIT1,HP-27,TRANSCRIPTIONFACTORSP1,GLYCERALDEHYDE-3-PHOSPHATEDEHYDROGENASE,HEATSHOCKFACTOR1,HEATSHOCK70KDAPROTEIN1A,SERUMALBUMINPREPROTEIN,SERUMALBUMINPREPROPROTEIN
	Ending time - 2021-06-07 06:44:23 - Neotamias

The total data results are also in the the same file folder with the A_Summary.txt in the file A_Total_Table.dat. This file has the following format. 

uniqueID	DB	ID	Accession	Genus_Species	Genus	Species	BIN_OTU	Gene	Sequence	Flags
2	GenBank	MT262677	MT262677	Tamias minimus	Tamias	minimus	GenBank	CYTOCHROMEB	AGCT	Wrong_Taxa

The column 'Flags' holds several different potential variables including, no_marker, no_taxa, no_seq, name_issue, taxa_digits, taxa_punct, wrong_taxa, or '-'. The '-' result indicates that the record was suitable for inclusion in the data set with no potential flags. 

## Using the create_fastas() function to select the records of interest to create a fasta file of molecular sequence data.
Using the A_summary.txt file from the auto_seq_download() function a table of target genera and associated names for target molecular markers needs to be manually assembled. There will most likely be multiple names for the same molecular marker due to variation in naming conventions. The format to construct the file is shown below. The file used in this example is included in the Example file folder and is named Chipmunk_marker_table.dat. 

Eutamias	Eutamias	Tamias	Tamias	Neotamias	Neotamias
CYTOCHROMEB	CYTOCHROMECOXIDASESUBUNITI	CYTB	COI-5P	CYTOCHROMEB	CYTOCHROMECOXIDASESUBUNITI
CYTOCHROME-B	CYTOCHROMEOXIDASESUBUNIT1	CYTOCHROMEB	CYTOCHROMECOXIDASESUBUNIT1	CYTOCHROME-B	CYTOCHROMEOXIDASESUBUNIT1
CYTOCHROMEBOXIDASE	CYTOCHROMEOXIDASESUBUNITI	CYTOCHROME-B	CYTOCHROMECOXIDASESUBUNITI	CYTOCHROMEBOXIDASE	CYTOCHROMEOXIDASESUBUNITI
		CYTOCHROMEBOXIDASE	CYTOCHROMEOXIDASESUBUNIT1		
			CYTOCHROMEOXIDASESUBUNITI		

The A_Total_Table.dat and the constructed Chipmunk_marker_table.dat were selected when prompted after running the create_fastas() function. For this example when running the create_fastas() function no arguments were used as all of the default arguments were wanted. These defaults would not include any records that had a flag, indicating a potential non-target or problematic sequence. The flags in question are no_marker, no_taxa, no_seq, name_issue, taxa_digits, taxa_punct, wrong_taxa. Only records with '-' were used to construct the fasta files. The following is the output on the R screen for this example.

	> create_fastas()
	Please select the total tables file.  Hit enter key to continue...
	Please select the file with genus and the list of molecular markers of interest. Hit enter key to continue...
	"Eutamias - Unable to create fasta file, no records for this genus and molecular markers"
	"Eutamias - Unable to create fasta file, no records for this genus and molecular markers"
	"Complete, please see the file location with your input table for results."
	"Complete, please see the file location with your input table for results."
	"Neotamias - Unable to create fasta file, no records for this genus and molecular markers"
	"Neotamias - Unable to create fasta file, no records for this genus and molecular markers"

The output from the create_fastas() function is the creation of fasta files. In this example two fasta files were created Tamias_COI-5P.fas and Tamias_CYTB.fas and they are both included in the 'Example' file folder. After creating the fasta files of interest using the create_fastas() function and the above steps, they need to be manually sorted and placed in to filefolders for each target molecular marker. In the example they were placed in to the COI and CytB. 

## Align the fasta files to a reference sequence. 
For this example, we will align the COI sequences to the closely related taxa Sciurus carolinensis. The reference fasta file used in this example is included in the 'Example' file folder. Please note that the header of the reference sequence needs to reflect the same format as those in the fasta files created in the create_fastas() function. Also, there can be no spaces in the naming for the reference fasta file. This format is shown below.

	>ABMC288-05||Sciurus|carolinensis|JF457099|COI-5P

Once running the align_to_ref() the file folder with the target sequences and the reference sequence will be prompted to be selected. In addition, the external to R program MAFFT and the file folder location of the program on the local computer will be requested. There are two arguments that can be include when calling the align_to_ref() function. The first is the 'pigl' argument. This indicates the percent of records with a gap at a particular position allowed before removing the records responsible for creating that gap. For this example the default of 0.95% was used. The second variable is the op variable. This indicate the the opening gap penalty used int he MAFFT alignment. For conserved regions this could be set to a higher value. The default of 1.53 was used for this example. However, since this example is using the COI-5P region, a gene which codes to a protein, this could have been increased to something like 10. The following is the output on the R screen for this example.

	> align_to_ref()
	Choose the folder location where your fasta files to be aligned are located. Hit enter key to continue...
	Choose your fasta reference file (note this must be a trimmed file with all sequences of the same length and no leading or trailing gap characters). Hit enter key to continue...
	Choose the folder location where the MAFFT (.bat file) is located. Hit enter key to continue...
	"C:\\A_MACER\\Seq_auto_dl_063034_Jun_07\\Total_Tables\\COI/Tamias_COI-5P.fas at 2021-06-17 05:47:14"
	"Length of the fasta file reported in nucleotides 1539 and the number of records 217"
	"Length of the aligned and trimmed multiple sequence alignment (MSA) in nucleotides 600 and the number of records 189"
	"Output is located in the target directory in the subfolders MAFFT and MAFFT_trimmed"

The output from the align_to_ref() function will be located in the folder with the target fasta files that was selected by the user during the function run. There will be three items in this folder. 
  1. a file folder with all MAFFT alignments to the reference sequence in MSA fasta files for each input target fasta file.
  2. a file folder with the aligned files from the first folder but trimmed to the reference sequence length and not including the reference sequence.
  3. a file with the log of the MAFFT alignment

## Evaluate the downloaded sequence for quality.

The final step in the process is to use the aligned and trimmed sequences and check the records for suitable inclusion in the dataset. The barcode_clean() function will be utilized to evaluate the resulting multiple sequence alignments (MSA). There are two arguments for the barcode_clean() function. The AA_code argument accepts a numerical value which identifies the amino acid coding matrix to be used. If this function is set to 0 then no check for amino acid coding will be completed. If this is used the input MSA must be in reading frame to obtain accurate results. If the reference sequenced used to construct the MSA was in reading frame then this should not be an issue. We utilized this argument in this example as the default value is 5, or invertebrate mitochondrial translation matrix, where we require 2 the vertebrate mitochondrial matrix. The second argument is the AGCT_only argument. This will eliminate sequences with characters other than AGCT (in other words remove all characters with uncertain IUPAC nucleotide coding). This argument is set to 1, on or in use. We did not use this argument as the default of on was what was required. When prompted the MAFFT_trimmed file folder generated by the align_to_ref() function was selected. The following is the output on the R screen for this example.

  > barcode_clean(AA_code=2)
  Choose the folder location where your input files are located. Hit enter key to continue...
  [1] "Start time...  2021-06-17 06:24:49"
            Species Initial Dereplicate AGCT AA Genus_Outlier Species_Outlier       Intraspecific       Interspecific Barcode_Gap
  1         amoenus      21          12   12 12            12              12  0.0666666666666667 0.00166666666666667          NO
  2         canipes       6           6    6  6             6               6 0.00833333333333333               0.025         YES
  3   cinereicollis      10          10   10 10            10              10  0.0116666666666667                   0          NO
  4        dorsalis      12          12   12 12            12              12  0.0216666666666667                   0          NO
  5         minimus       2           2    2  2             2               2                0.14 0.00166666666666667          NO
  6  quadrivittatus      12          12   12 12            12              12  0.0233333333333333                   0          NO
  7      ruficaudus       2           2    2  2             2               2  0.0333333333333333 0.00166666666666667          NO
  8           rufus       6           6    6  6             6               6                   0  0.0116666666666667         YES
  9       sibiricus      71          51   51 51            51              49  0.0216666666666667               0.115         YES
  10       striatus      36          19   19 19            19              19  0.0216666666666667 0.00166666666666667          NO
  11       umbrinus      11          11   11 11            11              11  0.0216666666666667                   0          NO
  [1] "Start time...  2021-06-17 06:24:49  and end time...  2021-06-17 06:24:50"

The results of this function will be placed in the MAFFT_trimmed file folder. The results are in three different files for each input fasta file and a single log file for the running of the function. The log file will be named with the format A_Clean_File_YYY-MM-DDtttttt and the file generated by this example is A_Clean_File_2021-06-17062449.dat. In addition to this file, the three files for each input fasta are:
  1. a data table file indicated with '_data_table.dat' appended to the end of the input file name for designation
  2. a distance matrix file with the matrix used in the function, indicated with '_dist_matrix.dat' appended to the end of the input file name for designation
  3. a fasta file with all outliers removed, in a file named using the input file name appended with '_no_outliers.fas'
The fasta file generated using this function will remove all records identified as non-unique or potentially inaccurate. The files are first dereplicated and duplicate records with the same GenBank accession are removed. Then like the auto_seq_download() function, the barcode_clean() function flags the remaining records. The column 'Flags' holds several different potential variables including, non_AGCT, Stop_Codon, Genus_Outlier, Species_Outlier, and '-'. The '-' result indicates that the record was suitable for inclusion in the data set with no potential flags. The generated fasta file removes all records with flags other than '-'. If the outcoming fasta files are too strict in their removal of records, for certain research purposes, the results of the function and the flags associated with each of the records can be viewed in the data table file. Finally, the inclusion of the distance matrix used to calculate the outliers for the function is also included as this file may want to be used for other analyses outside the scope of MACER.
>>>>>>> 9fb0bd654917b6bb749a7331c8260020217d5e32
