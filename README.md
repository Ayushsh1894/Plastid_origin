**Scirpt one: "Tree_sorting.txt"**

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
