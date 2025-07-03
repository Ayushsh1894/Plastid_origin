**Dependencies**

bash
# Python packages

**pip install biopython ete3**

# R packages  

**install.packages("PhySortR")**

#################################################################


**Workflow Execution**

Step 1: **Tree Sorting**

bash# Run in RStudio or R environment

Rscript Tree_sorting.txt

Step 2: **Sequence Extraction**
# **For Apicomplexa/Chromerida**
python extract_fasta_Apicomplexa.py /path/to/Chromeride.faa

# For Dinoflagellates  
python extract_fasta_Dinoflagellata.py /path/to/Dinoflagellate.faa

Step 3: **Transfer Analysis**

# **Chromerida transfer detection**
python Transfers_analysis_Apicomplexa.py

# Dinoflagellate transfer detection  
python Transfer_analysis_Dinoflagellata.py


#####################################################

Detailed explaination


**Script 1: "Tree_sorting.txt"**

Purpose: This R script performs phylogenetic tree sorting using the PhySortR package to analyze and categorize trees based on specific taxonomic groups.
Functionality:

Library: Uses the **PhySortR** package for phylogenetic tree analysis
Target Groups: Focuses on three taxonomic groups - Chromerida, Cryptophyta, and Cryptophyta (appears to be a duplicate entry)
Filtering: Sets a minimum proportion threshold of 0.1 (10%) for target groups to be included in the analysis
Input/Output:

Reads tree files from the Api_crypto_tree_files directory
Outputs sorted results to the Chrom_Crypto_crypto directory


Processing Mode: Uses mode "c" (likely for classification or clustering)
Clade Sorting: Configured to sort clades with "E" parameter
File Format: Processes files with .alg.fas.treefile extension (aligned FASTA tree files)

Taxonomic Context: The script includes comments listing various taxonomic groups of interest:

Stramenopiles
Apicomplexa
Chromerida (Vbra and Cvel variants)
Dinoflagellates
Haptophyta
Cryptophyta
Rhodophyta

##################################################


**Script 2: **extract_fasta_Apicomplexa.py****

Purpose: This script performs phylogenetic analysis and sequence extraction for Apicomplexa and Chromerida lineages from tree files and FASTA sequences.
Key Functionality:

#Tree Processing:

Loads phylogenetic trees using the ete3 library
Identifies subtrees containing both Chromerida and target taxonomic groups
Validates subtree structure to ensure proper sister-group relationships
Handles special cases like Archaea_ASGARD as a composite target group

#Sequence Extraction:

Parses FASTA files to extract protein sequences
Maps tree leaf names to sequence headers using pattern matching
Extracts sequences based on orthologous group (OG) identifiers
Handles complex ID parsing from headers like >Chromerida-Cvel_XXXXX_tX_bY_OGZ

#Data Processing Pipeline:

Subtree Identification: Finds valid sister-group relationships between Chromerida and target taxa
Sequence Mapping: Links phylogenetic tree leaves to corresponding protein sequences
Duplicate Removal: Filters sequences to keep only one per branch-OG combination
Organism Separation: Splits sequences by organism type (Cvel vs Vbra)

#Output Structure:

merged_fasta/: Combined sequences by taxonomic group
filtered_fasta/: Deduplicated sequences
split_fasta_files/: Organism-specific sequence files
Individual subtree files (.nw format) and sequence ID lists

Input Requirements:

Directory structure: Apicomplexa-[TARGET_GROUP]/ containing .treefile files
FASTA file with Chromerida protein sequences
Can accept input via command line argument or interactive prompt

#################################################

**Script 3: **extract_fasta_Dinoflagellata.py****

Purpose: This script performs phylogenetic analysis and sequence extraction for Dinoflagellate lineages from tree files and FASTA sequences, with enhanced error handling and progress monitoring.
Key Functionality:
Tree Processing:

Loads phylogenetic trees using the ete3 library with timeout protection (60 seconds)
Identifies subtrees containing both Dinoflagellates and target taxonomic groups
Validates subtree structure to ensure proper sister-group relationships
Handles special cases like Archaea_ASGARD as a composite target group ["Archaea", "ASGARD"]
Implements robust error handling with comprehensive exception management

#Sequence Extraction:

Parses FASTA files to extract protein sequences with error recovery
Maps tree leaf names to sequence headers using pattern matching
Extracts sequences based on orthologous group (OG) identifiers
Handles complex ID parsing from headers like >Dinoflagellate-Species_ID_bX_OGY
Provides progress tracking for large datasets (reports every 10 files)

#Data Processing Pipeline:

Subtree Identification: Finds valid sister-group relationships between Dinoflagellates and target taxa
Sequence Mapping: Links phylogenetic tree leaves to corresponding protein sequences
Duplicate Removal: Filters sequences to keep only one per branch-OG combination
Quality Control: Validates parsing success and reports statistics
Progress Monitoring: Creates detailed execution logs and timing reports

#Output Structure:

merged_dino_fasta/: Combined sequences by taxonomic group
filtered_dino_fasta/: Deduplicated sequences
progress_tracking/: Individual progress files for monitoring
Individual subtree files (.nw format) and sequence ID lists
dino_extraction_summary.txt: Comprehensive execution report with statistics

Enhanced Features:

Timeout Protection: Prevents hanging on computationally intensive large trees
Real-time Progress: Shows current processing status and estimated completion
Error Recovery: Continues processing even when individual files fail
Detailed Logging: Comprehensive summary with timing, success rates, and error documentation
Memory Optimization: Efficient processing for large-scale datasets

Input Requirements:

Directory structure: Dinoflagellates-[TARGET_GROUP]/ containing .treefile files
FASTA file with Dinoflagellate protein sequences
Can accept input via command line argument, interactive prompt, or stdin (for batch jobs)
Automatically discovers and processes all relevant directories

Performance Features:

Batch Processing: Handles multiple target groups simultaneously
Progress Files: Individual tracking for each directory being processed
Statistics Reporting: Detailed metrics on processing success and failure rates
Time Tracking: Monitors execution time for optimization and planning

##################################################################

**Script 4: Transfers_analysis_Apicomplexa.py**

Purpose: Parallel analysis detecting gene transfer events involving Chromerida (closely related to Apicomplexa).
Core Functionality: Nearly identical to the Dinoflagellates script with key adaptations:
Target Group Substitution:

Analyzes Chromerida instead of Dinoflagellates
Processes Apicomplexa-* directories instead of Dinoflagellates-*
Uses Chromeride.faa as the reference sequence database

Analysis Patterns:

Type A: (Chromerida + Target Group) sister to (Pure Target Group)
Type B: (Chromerida + Target Group) sister to (Pure Chromerida)

Key Functions (Chromerida-specific versions):

find_chromerida_targetgroup_sister_to_pure_chromerida(): Detects Chromerida-specific transfer patterns
extract_chromerids_ids(): Extracts Chromerida sequence identifiers
Same tree topology validation functions with Chromerida prefixes

Input/Output Structure:

Input: Apicomplexa-[TARGET_GROUP]/ directories with .treefile files
Reference: /home/users/ayush/Ayush/New/DBs/Chromeride.faa
Output: Separate_Chromerida_Analysis_Results/ with same substructure as Dinoflagellates script

Common Features of Both Scripts:
Comprehensive Target Group Coverage:

Analyzes 19 different taxonomic groups including major eukaryotic and prokaryotic lineages
Handles complex relationships like Archaea/ASGARD combinations

Robust Error Handling:

Dependency checking for Biopython and ete3
File validation and graceful error recovery
Detailed logging and progress reporting

Statistical Analysis:

Counts valid trees for each relationship type
Tracks unique sequences involved in transfers
Generates summary statistics and reports
Prevents double-counting in merged analyses

File Management:

Automated directory creation and organization
Tree file copying to appropriate result directories
FASTA sequence extraction and writing
Comprehensive summary report generation

Research Application:
These scripts are designed for studying horizontal gene transfer in protist evolution, particularly focusing on:

Chromerida-Apicomplexa relationships and their connections to other lineages
Dinoflagellate evolutionary history and potential gene acquisitions
Comparative analysis of transfer patterns between these two related groups

############################

**Script 5: Transfer_analysis_Dinoflagellata.py**

Purpose: This script analyzes phylogenetic trees to identify potential gene transfers events involving Dinoflagellates by detecting specific sister-group relationships that suggest gene transfer.
Core Analysis Types:

1. Sister-Group Relationship Detection:

Type A: (Dinoflagellates + Target Group) sister to (Pure Target Group)
Type B: (Dinoflagellates + Target Group) sister to (Pure Dinoflagellates)

Key Functions:
Tree Topology Analysis:

is_pure_clade(): Validates that a clade contains only sequences from one taxonomic group
has_both_groups(): Confirms presence of sequences from two specified groups
find_mixed_pure_sister_clades(): Identifies sister relationships between mixed and pure clades
find_dinoflagellates_targetgroup_sister_to_pure_dinoflagellates(): Finds specific Dinoflagellates transfer patterns

Special Handling:

Archaea/ASGARD Processing: Treats Archaea and ASGARD as separate groups but creates merged results for comprehensive analysis
Composite Group Support: Handles Archaea_ASGARD as a combined target group

Processing Pipeline:

Directory Processing: Analyzes all Dinoflagellates-* folders containing tree files
Tree Validation: Checks each .treefile for valid sister-group relationships
Sequence Extraction: Extracts Dinoflagellates sequences involved in potential transfers
Results Organization: Creates separate output directories for each relationship type
Summary Generation: Produces comprehensive reports with statistics

Output Structure:

Valid_Dinoflagellates_Sister_Trees/: Trees showing Dinoflagellates transfer patterns
Valid_[TARGET]_Sister_Trees/: Trees for each target group analysis
Valid_Archaea_ASGARD_Merged_Sister_Trees/: Combined Archaea/ASGARD results
Extracted_Dinoflagellates_[TARGET]_sequences.faa: Sequence files for each relationship
Sister_Relationships_Summary.txt: Detailed analysis report

##############################################






