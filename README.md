**script 1: "Tree_sorting.txt"**

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

#################################################################################################
**Script 2: **extract_fasta_Apicomplexa.py****
Purpose: This script performs phylogenetic analysis and sequence extraction for Apicomplexa and Chromerida lineages from tree files and FASTA sequences.
Key Functionality:
Tree Processing:

Loads phylogenetic trees using the ete3 library
Identifies subtrees containing both Chromerida and target taxonomic groups
Validates subtree structure to ensure proper sister-group relationships
Handles special cases like Archaea_ASGARD as a composite target group

Sequence Extraction:

Parses FASTA files to extract protein sequences
Maps tree leaf names to sequence headers using pattern matching
Extracts sequences based on orthologous group (OG) identifiers
Handles complex ID parsing from headers like >Chromerida-Cvel_XXXXX_tX_bY_OGZ

Data Processing Pipeline:

Subtree Identification: Finds valid sister-group relationships between Chromerida and target taxa
Sequence Mapping: Links phylogenetic tree leaves to corresponding protein sequences
Duplicate Removal: Filters sequences to keep only one per branch-OG combination
Organism Separation: Splits sequences by organism type (Cvel vs Vbra)

Output Structure:

merged_fasta/: Combined sequences by taxonomic group
filtered_fasta/: Deduplicated sequences
split_fasta_files/: Organism-specific sequence files
Individual subtree files (.nw format) and sequence ID lists

Input Requirements:

Directory structure: Apicomplexa-[TARGET_GROUP]/ containing .treefile files
FASTA file with Chromerida protein sequences
Can accept input via command line argument or interactive prompt

###################################################################################################
**Script 3: **extract_fasta_Dinoflagellata.py****
Purpose: Similar to the Apicomplexa script but specialized for analyzing Dinoflagellate phylogenetic relationships with other taxonomic groups.
Key Differences from Apicomplexa Script:
Enhanced Error Handling:

Implements timeout protection (60 seconds) for large tree processing
Comprehensive exception handling with detailed error logging
Progress tracking with intermediate status reports

Performance Optimizations:

Batch processing with progress indicators
Memory-efficient sequence handling
Timeout mechanisms to prevent hanging on complex trees

Robust Processing:

Creates detailed summary files with execution statistics
Progress tracking files for monitoring long-running jobs
Comprehensive logging of processing times and results

Target Focus:

Searches for subtrees containing both Dinoflagellates and target groups
Uses Dinoflagellate- prefix instead of Chromerida- for sequence identification
Handles the same special cases (e.g., Archaea_ASGARD)

Output Structure:

merged_dino_fasta/: Combined Dinoflagellate sequences by group
filtered_dino_fasta/: Deduplicated sequences
progress_tracking/: Individual progress files for each directory
dino_extraction_summary.txt: Comprehensive execution summary

Enhanced Monitoring:

Real-time progress reporting during execution
Detailed timing information for each processing step
Summary statistics for troubleshooting and optimization

Common Features of Both Scripts:

Process multiple target taxonomic groups simultaneously
Handle complex phylogenetic tree structures
Extract and organize sequences based on orthologous groups
Remove duplicate sequences while preserving taxonomic information
Support both single target groups and composite groups (lists)
Generate standardized output directory structures
#####################################################################################################


