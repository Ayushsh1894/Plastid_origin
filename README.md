# Plastid_origin
Independent origins of apicomplexan and dinoflagellate plastids, with contributions from multiple sources

####################################################
**"Tree_sorting.txt"**



This R script uses the PhySortR package to automatically sort and filter phylogenetic trees based on specific taxonomic groups of interest.
Purpose
The script identifies and extracts phylogenetic trees that contain representatives from target taxonomic groups (Chromerida and Cryptophyta) from a larger collection of tree files.
Key Parameters

Target Groups: Chromerida, Cryptophyta, Cryptophyta (specified taxonomic clades)
Minimum Proportion: 0.1 (trees must contain at least 10% of sequences from target groups)
Input Directory: Api_crypto_tree_files/ (contains original tree files)
Output Directory: Chrom_Crypto_crypto/ (filtered trees are saved here)
Mode: "c" (conservative filtering mode)
Clades Sorted: "E" (sorts by evolutionary relationships)
File Extension: .alg.fas.treefile (specific tree file format)

Taxonomic Groups Referenced
The script includes comments listing various protist groups that appear to be part of a larger phylogenetic study:

Stramenopiles
Apicomplexa
Chromerida (Vbra and Cvel strains)
Dinoflagellates
Haptophyta
Cryptophyta
Rhodophyta


Usage
This script is designed to run in RStudio and requires the PhySortR package to be installed. It's particularly useful for large-scale phylogenetic analyses where you need to quickly identify trees containing specific taxonomic groups of interest.

####################################################

**"extract_fasta_seq.py**"


This Python script performs comprehensive phylogenetic analysis to identify and extract potential gene transfer events between Chromerida and other taxonomic groups.
Purpose
The script processes phylogenetic trees to identify subtrees containing both Chromerida sequences and sequences from target taxonomic groups, then extracts and processes the corresponding protein sequences for further analysis.
Key Features
1. Tree Analysis

Loads phylogenetic trees using the ETE3 toolkit
Identifies subtrees containing both Chromerida and target group sequences
Validates subtree structure to ensure meaningful evolutionary relationships
Extracts and saves relevant subtrees in Newick format

2. Sequence Processing

Reads FASTA files containing protein sequences
Extracts sequences corresponding to identified phylogenetic relationships
Filters sequences to remove duplicates based on branch-OG combinations
Splits sequences by organism (Cvel and Vbra strains)

3. Data Organization

Creates organized directory structure for outputs:

merged_fasta/: Combined sequences by taxonomic group
filtered_fasta/: Deduplicated sequences
split_fasta_files/: Organism-specific sequence files



Dependencies

ete3: For phylogenetic tree manipulation
os, re, sys: Standard Python libraries for file operations
collections.defaultdict: For efficient data grouping

Input Requirements

Directory structure with Apicomplexa-[GROUP] folders containing .treefile files
FASTA file with Chromerida protein sequences
Tree files should follow the naming convention including OG (Orthologous Group) numbers

Target Groups Supported

Standard taxonomic groups (e.g., Cryptophyta, Stramenopiles)
Special handling for complex groups like Archaea_ASGARD
Automatic detection of target groups from directory names

Output Files

Subtree files: .nw format containing relevant phylogenetic relationships
Sequence ID files: Text files listing sequences in each subtree
Merged FASTA files: Combined sequences by taxonomic group
Filtered FASTA files: Deduplicated sequences
Split FASTA files: Separate files for Cvel and Vbra organisms

Usage
bashpython Transfers_analysis.py [path_to_fasta_file]
The script can accept the FASTA file path as a command-line argument or through interactive input.
Applications
This script is particularly useful for:

Identifying potential horizontal gene transfer events
Comparative genomics studies
Phylogenetic analysis of specific gene families
Evolutionary studies of protist lineages



######################################################################
**"Transfers_analysis.py"**

This Python script performs advanced phylogenetic analysis to identify specific sister-group relationships and extract corresponding protein sequences, particularly focusing on potential horizontal gene transfer events involving Chromerida.
Purpose
The script analyzes phylogenetic trees to identify two key evolutionary patterns:

Mixed clades containing both Chromerida and target group sequences that are sister to pure target group clades
Chromerida relationships where mixed (Chromerida + target group) clades are sister to pure Chromerida clades

Key Features
1. Sister Relationship Detection

Pure Clade Validation: Ensures clades contain only sequences from a single taxonomic group
Mixed Clade Analysis: Identifies clades containing exactly two taxonomic groups (no contamination)
Sister Group Identification: Finds specific evolutionary relationships between mixed and pure clades

2. Special Group Handling

Archaea_ASGARD: Special logic for complex archaeal groups (Archaea + ASGARD)
Merged Analysis: Combines results from related groups to avoid double-counting
Flexible Target Groups: Supports 19 different taxonomic groups

3. Sequence Extraction

ID Matching: Extracts Chromerida sequences corresponding to identified phylogenetic patterns
FASTA Processing: Uses BioPython for robust sequence file handling
Organized Output: Creates separate FASTA files for each target group

4. Data Organization

Folder-Specific Processing: Handles multiple Apicomplexa-[GROUP] input directories
Structured Output: Creates organized directory hierarchy for results
Valid Tree Copying: Copies only trees with valid relationships to output folders

Dependencies

BioPython: For FASTA file parsing and sequence manipulation
ETE3: For phylogenetic tree analysis and traversal
Standard Libraries: os, sys, datetime, shutil, collections

Target Groups Analyzed
The script analyzes 19 taxonomic groups:

Dinoflagellates, Rhizaria, Stramenopiles
Glaucocystophyceae, Cryptophyta, Haptophyta
Rhodophyta, Discoba, Rhodelphis
Opisthokonta, Metamonada, Evosea
Archaea, ASGARD, Archaea_ASGARD
Proteobacteria, Viridiplantae, Cyanobacteria
Chromerida

Input Requirements

Directory Structure: Multiple Apicomplexa-[GROUP] folders containing .treefile files
FASTA Database: Chromerida protein sequence database (Chromeride.faa)
Tree Format: Newick format trees with proper taxonomic prefixes

Output Structure
Separate_Chromerida_Analysis_Results/
├── Apicomplexa-[GROUP1]/
│   ├── Valid_[GROUP]_Sister_Trees/
│   ├── Valid_Chromerida_Sister_Trees/
│   ├── Extracted_Chromerida_[GROUP]_sequences.faa
│   └── Sister_Relationships_Summary.txt
└── Apicomplexa-[GROUP2]/
    └── ...
Key Functions
1. Tree Analysis Functions

is_pure_clade(): Validates taxonomic purity of clades
has_both_groups(): Checks for presence of both target groups
find_mixed_pure_sister_clades(): Identifies main sister relationships
find_chromerida_targetgroup_sister_to_pure_chromerida(): Finds Chromerida-specific patterns

2. Processing Functions

process_single_folder(): Handles individual input directories
extract_chromerids_ids(): Extracts Chromerida sequence identifiers
summary_report(): Generates comprehensive analysis reports

Applications
This script is particularly valuable for:

Horizontal Gene Transfer Studies: Identifying potential HGT events
Evolutionary Analysis: Understanding sister-group relationships
Comparative Genomics: Analyzing gene family evolution
Phylogenetic Validation: Confirming evolutionary hypotheses

Usage
bashpython extract_fasta_seq.py
The script automatically processes all Apicomplexa-* directories in the current working directory and requires the Chromerida FASTA database to be available at the specified path.
Output Reports
Each processed folder generates a detailed summary report including:

Number of valid trees found per target group
Unique Chromerida sequences identified
Overall statistics and relationships
Processing timestamps and metadata


################################################################################################################
**COMPLETE EXPLAINATION**

# Plastid_origin

A comprehensive bioinformatics pipeline for analyzing horizontal gene transfer (HGT) events and evolutionary relationships between Chromerida and other taxonomic groups, with particular focus on plastid origin and evolution.

## Overview

This repository contains a complete workflow for phylogenetic analysis of potential horizontal gene transfer events involving Chromerida, a key group in understanding plastid evolution. The pipeline processes large-scale phylogenetic datasets to identify, extract, and analyze evolutionary relationships that may indicate gene transfer events.

## Pipeline Workflow

The analysis follows a systematic three-step approach:

```
1. Tree Sorting (R) → 2. Sister Relationship Analysis (Python) → 3. Transfer Analysis (Python)
```

### Step 1: Initial Tree Filtering
**Script**: `Tree_sorting.txt` (R/RStudio)
- Filters phylogenetic trees containing target taxonomic groups
- Uses PhyloSortR package for automated tree selection
- Identifies trees with sufficient representation of Chromerida and target groups
- **Input**: Directory of phylogenetic tree files (`.treefile` format)
- **Output**: Filtered trees in organized directories

### Step 2: Sister Relationship Analysis & Sequence Extraction
**Script**: `extract_fasta_seq.py` (Python)
- Identifies specific sister-group relationships in phylogenetic trees
- Detects patterns: (Chromerida + Target Group) sister to (Pure Target Group)
- Extracts corresponding protein sequences for further analysis
- **Input**: Filtered trees from Step 1, Chromerida FASTA database
- **Output**: Valid trees, extracted sequences, relationship summaries

### Step 3: Comprehensive Transfer Analysis
**Script**: `Transfers_analysis.py` (Python)
- Performs detailed analysis of potential HGT events
- Processes subtrees containing both Chromerida and target sequences
- Organizes sequences by organism and removes duplicates
- **Input**: Processed trees and sequences from Step 2
- **Output**: Organized sequence files, filtered datasets, split organism files

## Key Features

### Multi-Scale Analysis
- **Phylogenetic Scale**: Tree topology analysis for evolutionary relationships
- **Sequence Scale**: Protein sequence extraction and organization
- **Taxonomic Scale**: Analysis across 19+ taxonomic groups

### Robust Data Processing
- **Quality Control**: Multiple validation steps for tree and sequence data
- **Error Handling**: Comprehensive error checking and logging
- **Scalability**: Designed for large-scale phylogenomic datasets

### Specialized Algorithms
- **Sister Relationship Detection**: Advanced algorithms for identifying specific evolutionary patterns
- **Contamination Detection**: Ensures taxonomic purity of analyzed clades
- **Duplicate Removal**: Sophisticated filtering based on orthologous group assignments

## Target Taxonomic Groups

The pipeline analyzes evolutionary relationships between Chromerida and:

**Primary Eukaryotic Groups:**
- Dinoflagellates, Stramenopiles, Cryptophyta
- Haptophyta, Rhodophyta, Glaucocystophyceae
- Rhizaria, Discoba, Rhodelphis
- Opisthokonta, Metamonada, Evosea
- Viridiplantae

**Prokaryotic Groups:**
- Archaea, ASGARD, Proteobacteria, Cyanobacteria

**Special Categories:**
- Archaea_ASGARD (merged analysis)
- Chromerida (self-relationships)

## Technical Requirements

### Software Dependencies
- **R/RStudio**: PhyloSortR package
- **Python 3.x**: BioPython, ETE3
- **System**: Unix/Linux environment (recommended)

### Input Data Structure
```
project_directory/
├── Api_crypto_tree_files/          # Original tree files
├── Apicomplexa-[GROUP1]/           # Group-specific directories
├── Apicomplexa-[GROUP2]/           # (created by pipeline)
├── Chromeride.faa                  # Chromerida protein database
└── analysis_scripts/               # This repository
```

### Output Data Structure
```
project_directory/
├── Chrom_Crypto_crypto/                    # Step 1 output
├── Separate_Chromerida_Analysis_Results/   # Step 2 output
│   ├── Apicomplexa-[GROUP]/
│   │   ├── Valid_[GROUP]_Sister_Trees/
│   │   ├── Extracted_Chromerida_sequences.faa
│   │   └── Sister_Relationships_Summary.txt
├── merged_fasta/                           # Step 3 output
├── filtered_fasta/
└── split_fasta_files/
```

## Usage Instructions

### Step 1: Tree Sorting
```r
# Run in RStudio
source("Tree_sorting.txt")
```

### Step 2: Sister Relationship Analysis
```bash
python extract_fasta_seq.py
```

### Step 3: Transfer Analysis
```bash
python Transfers_analysis.py [path_to_chromerida_fasta]
```

## Scientific Applications

### Horizontal Gene Transfer Studies
- Identification of potential HGT events between Chromerida and other groups
- Quantification of gene transfer frequencies across taxonomic boundaries
- Analysis of evolutionary patterns in plastid-related genes

### Phylogenetic Analysis
- Sister-group relationship validation
- Evolutionary tree topology analysis
- Taxonomic classification refinement

### Comparative Genomics
- Cross-taxonomic gene family analysis
- Orthologous group evolution tracking
- Protein sequence conservation studies

## Output Interpretation

### Tree Analysis Results
- **Valid Trees**: Trees containing meaningful evolutionary relationships
- **Sequence IDs**: Identifiers for sequences involved in potential HGT events
- **Subtree Files**: Extracted phylogenetic relationships in Newick format

### Sequence Analysis Results
- **Merged Sequences**: Combined datasets by taxonomic group
- **Filtered Sequences**: Deduplicated, high-quality sequence sets
- **Split Files**: Organism-specific sequence collections (Cvel, Vbra strains)

### Summary Reports
- **Relationship Statistics**: Quantitative analysis of sister-group patterns
- **Sequence Counts**: Number of sequences per taxonomic group
- **Quality Metrics**: Processing statistics and validation results

## Performance Considerations

### Memory Requirements
- **Tree Files**: ~100MB-1GB per analysis depending on dataset size
- **Sequence Files**: ~10MB-100MB for extracted sequences
- **RAM**: Minimum 8GB recommended for large datasets

### Processing Time
- **Step 1**: Minutes to hours depending on tree collection size
- **Step 2**: Hours for comprehensive multi-group analysis
- **Step 3**: Minutes to hours for sequence processing

## Troubleshooting

### Common Issues
1. **Missing Dependencies**: Ensure all required packages are installed
2. **File Path Errors**: Verify correct directory structure and file locations
3. **Memory Issues**: Consider processing smaller subsets for very large datasets
4. **Format Compatibility**: Ensure tree files are in proper Newick format

### Error Handling
- All scripts include comprehensive error checking
- Log files provide detailed information about processing steps
- Validation routines ensure data integrity throughout the pipeline

## Citation

If you use this pipeline in your research, please cite the relevant papers and acknowledge the tools used:
- PhyloSortR package
- ETE3 toolkit
- BioPython library

## Contributing

This pipeline is designed for research use in evolutionary biology and phylogenetics. For questions, improvements, or bug reports, please refer to the individual script documentation or contact the repository maintainer.

## License

This project is available for academic and research use. Please ensure proper attribution when using or modifying the code.
