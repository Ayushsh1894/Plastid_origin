from ete3 import Tree
import os
import re
import sys
import time
import traceback
from collections import defaultdict

def is_target_leaf(leaf_name, target_group):
    # Extract the taxonomy from folder name
    if isinstance(target_group, list):
        return (leaf_name.startswith("Dinoflagellate") or
                any(leaf_name.startswith(group) for group in target_group))
    else:
        return (leaf_name.startswith("Dinoflagellate") or
                leaf_name.startswith(target_group))

def extract_dino_subtrees(tree_file, target_group, timeout=60):
    """Extract subtrees with a timeout to prevent hanging on large trees"""
    print(f"Processing file: {tree_file}")
    start_time = time.time()
    
    try:
        tree = Tree(tree_file)
        print(f"Tree loaded, number of leaves: {len(tree.get_leaves())}")
        
        def is_valid_subtree(node):
            # Check if processing is taking too long
            if time.time() - start_time > timeout:
                raise TimeoutError(f"Processing timeout after {timeout} seconds")
                
            children = node.get_children()
            dino_child = None
            target_child = None
            other_children = []

            for child in children:
                if all(leaf.name.startswith("Dinoflagellate") for leaf in child.iter_leaves()):
                    dino_child = child
                # Handle both string and list target groups
                elif isinstance(target_group, list):
                    if all(any(leaf.name.startswith(group) for group in target_group) for leaf in child.iter_leaves()):
                        target_child = child
                # If target_group is a string, check if it matches
                elif all(leaf.name.startswith(target_group) for leaf in child.iter_leaves()):
                    target_child = child
                else:
                    other_children.append(child)

            if dino_child and target_child and len(other_children) == 0:
                return True
            elif dino_child and target_child and all(is_valid_subtree(child) for child in other_children):
                return True
            return False

        subtrees = []
        for node in tree.traverse(strategy="postorder"):
            if is_valid_subtree(node):
                subtrees.append(node)
                print(f"Found a suitable subtree with {len(node.get_leaves())} leaves")

        print(f"Total suitable subtrees found: {len(subtrees)}")
        return subtrees
        
    except TimeoutError as e:
        print(f"WARNING: {str(e)} for file {tree_file}. Skipping to next file.")
        return []
    except Exception as e:
        print(f"ERROR: Failed to process {tree_file}: {str(e)}")
        traceback.print_exc()
        return []

def read_fasta(filename):
    """Read a FASTA file with proper error handling"""
    sequences = {}
    try:
        current_id = None
        current_seq = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]  # Remove '>' but keep the full header
                    current_seq = []
                else:
                    current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)
        print(f"Successfully read FASTA file: {filename} with {len(sequences)} sequences")
        return sequences
    except Exception as e:
        print(f"ERROR: Failed to read FASTA file {filename}: {str(e)}")
        traceback.print_exc()
        return {}

def extract_dino_sequences(sorted_trees_dir, fasta_file):
    """Extract dinoflagellate sequences from FASTA file based on sequence IDs"""
    dino_sequences = {}
    sequences = read_fasta(fasta_file)
    
    if not sequences:
        print("WARNING: No sequences loaded from FASTA file")
        return dino_sequences
        
    # Count sequence ID files
    id_files = [f for f in os.listdir(sorted_trees_dir) if f.endswith("_sequence_ids.txt")]
    print(f"Found {len(id_files)} sequence ID files in {sorted_trees_dir}")
    
    # Track progress
    processed_files = 0
    
    # Collect all sequence ID files
    for filename in id_files:
        processed_files += 1
        if processed_files % 10 == 0:
            print(f"Progress: Processed {processed_files}/{len(id_files)} ID files")
            
        try:
            with open(os.path.join(sorted_trees_dir, filename), 'r') as f:
                for line in f:
                    line = line.strip()
                    og_match = re.search(r'(OG\d+)$', line)
                    if og_match:
                        og_number = og_match.group(1)
                        # Extract the ID part after "Dinoflagellate-" and before "_b"
                        short_id = line.split('-')[-1].split('_b')[0]

                        match_found = False
                        for header, sequence in sequences.items():
                            if short_id in header:
                                dino_sequences[line] = (header, sequence)
                                match_found = True
                                break
                                
                        if not match_found:
                            print(f"WARNING: No match found for sequence ID: {short_id}")
        except Exception as e:
            print(f"ERROR: Failed to process ID file {filename}: {str(e)}")
            continue

    print(f"Total sequences extracted: {len(dino_sequences)}")
    return dino_sequences

def parse_id(id_string):
    # Extract parts from the FASTA header
    match = re.match(r'>Dinoflagellate-(.+)_b(\d+)_OG(\d+)', id_string)
    if match:
        gene_name = match.group(1)
        branch = match.group(2)
        og = match.group(3)
        return gene_name, branch, og
    return None, None, None

def sort_and_filter_sequences(input_file, output_file):
    """Sort and filter sequences to remove duplicates"""
    try:
        sequences = defaultdict(list)
        skipped_lines = 0
        total_sequences = 0
        
        with open(input_file, 'r') as f:
            current_id = None
            current_sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if current_id and current_sequence:
                        gene_name, branch, og = parse_id(current_id)
                        if gene_name is None or branch is None or og is None:
                            print(f"Warning: Skipping sequence with ID {current_id}")
                            skipped_lines += 1
                        else:
                            key = (branch, og)
                            sequences[key].append((current_id, ''.join(current_sequence)))
                            total_sequences += 1
                    # Reset for new sequence
                    current_id = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
                    
            # Process the last sequence
            if current_id and current_sequence:
                gene_name, branch, og = parse_id(current_id)
                if gene_name is None or branch is None or og is None:
                    print(f"Warning: Skipping sequence with ID {current_id}")
                    skipped_lines += 1
                else:
                    key = (branch, og)
                    sequences[key].append((current_id, ''.join(current_sequence)))
                    total_sequences += 1
                    
        # Write output, keeping only one sequence per branch-OG combination
        with open(output_file, 'w') as f:
            for key in sorted(sequences.keys()):
                # Keep only the first sequence for each unique branch-OG combination
                f.write(sequences[key][0][0] + '\n')
                f.write(sequences[key][0][1] + '\n')
                
        return len(sequences), total_sequences, skipped_lines
        
    except Exception as e:
        print(f"ERROR: Failed to sort and filter sequences in {input_file}: {str(e)}")
        traceback.print_exc()
        return 0, 0, 0

def main():
    # Print start time
    start_time = time.time()
    print(f"Script started at: {time.ctime()}")
    
    # Default paths
    base_dir = os.getcwd()  # Current directory (dino_overlaps)

    # Check if path is provided as command line argument
    if len(sys.argv) > 1:
        default_fasta_file = sys.argv[1]
    else:
        # Try to read from stdin (for sbatch heredoc)
        try:
            default_fasta_file = input().strip()
        except EOFError:
            default_fasta_file = input("Enter the path to your Dinoflagellate.faa file: ")

    print(f"Using FASTA file: {default_fasta_file}")

    # Create directories for merged and filtered FASTA files
    merged_fasta_dir = os.path.join(base_dir, "merged_dino_fasta")
    filtered_fasta_dir = os.path.join(base_dir, "filtered_dino_fasta")
    progress_dir = os.path.join(base_dir, "progress_tracking")

    os.makedirs(merged_fasta_dir, exist_ok=True)
    os.makedirs(filtered_fasta_dir, exist_ok=True)
    os.makedirs(progress_dir, exist_ok=True)

    # Create a summary file
    summary_file = os.path.join(base_dir, "dino_extraction_summary.txt")
    
    with open(summary_file, 'w') as summary:
        summary.write(f"Dinoflagellate Sequence Extraction Summary\n")
        summary.write(f"======================================\n\n")
        summary.write(f"Started at: {time.ctime()}\n")
        summary.write(f"FASTA file: {default_fasta_file}\n\n")
        
        # Create a dictionary to group sequences by their group
        group_sequences = {}
        
        # Get list of directories to process
        dirs_to_process = [d for d in os.listdir(base_dir) if d.startswith("Dinoflagellates-") and os.path.isdir(os.path.join(base_dir, d))]
        summary.write(f"Found {len(dirs_to_process)} directories to process:\n")
        for d in dirs_to_process:
            summary.write(f"  - {d}\n")
        summary.write("\n")

        # Process each directory (Dinoflagellates-*)
        for i, subdir in enumerate(dirs_to_process, 1):
            # Mark directory as in progress
            progress_file = os.path.join(progress_dir, f"{subdir}.progress")
            with open(progress_file, 'w') as pf:
                pf.write(f"Started processing at {time.ctime()}\n")
                
            current_dir = os.path.join(base_dir, subdir)
            dir_start_time = time.time()

            # Get the target group from directory name
            target_group = subdir.split("-")[1]

            # Special case for Archaea_ASGARD
            if target_group == "Archaea_ASGARD":
                target_group = ["Archaea", "ASGARD"]

            # For organizational purposes
            group_name = target_group.lower() if not isinstance(target_group, list) else "archaea_asgard"

            print(f"\nProcessing directory {i}/{len(dirs_to_process)}: {subdir}")
            print(f"Target group: {target_group}")
            
            summary.write(f"Directory {i}/{len(dirs_to_process)}: {subdir}\n")
            summary.write(f"  Target group: {target_group}\n")

            # Process each tree file in the directory
            tree_files = [f for f in os.listdir(current_dir) if f.endswith(".treefile")]
            tree_files_count = len(tree_files)
            
            if tree_files_count == 0:
                print(f"No .treefile files found in {current_dir}")
                summary.write(f"  No .treefile files found\n\n")
                continue
                
            print(f"Found {tree_files_count} tree files to process")
            summary.write(f"  Found {tree_files_count} tree files\n")
            
            subtrees_found = 0
            processed_trees = 0
            
            for j, filename in enumerate(tree_files, 1):
                print(f"\nProcessing file {j}/{tree_files_count}: {filename}")
                tree_file = os.path.join(current_dir, filename)
                
                try:
                    subtrees = extract_dino_subtrees(tree_file, target_group)
                    processed_trees += 1
                    
                    # Extract OG number from filename
                    og_match = re.search(r'(OG\d+)', filename)
                    og_number = og_match.group(1) if og_match else "UnknownOG"

                    if subtrees:
                        subtrees_found += len(subtrees)
                        print(f"Subtrees for {filename} containing Dinoflagellate and {target_group}:")
                        for i, subtree in enumerate(subtrees, 1):
                            sequence_ids = [leaf.name for leaf in subtree.iter_leaves() if is_target_leaf(leaf.name, target_group)]
                            print(f"Subtree {i}:")
                            print(sequence_ids)
                            print()

                            # Save the subtree
                            outfile = os.path.join(current_dir, f"{og_number}_dino_{str(target_group).lower()}_subtree_{i}.nw")
                            subtree.write(outfile=outfile, format=1)
                            print(f"Subtree saved to: {outfile}")

                            # Save sequence IDs to a file
                            id_file = os.path.join(current_dir, f"{og_number}_b{i}_dino_{str(target_group).lower()}_sequence_ids.txt")
                            with open(id_file, "w") as f:
                                for seq_id in sequence_ids:
                                    f.write(f"{seq_id}_b{i}_{og_number}\n")
                            print(f"Sequence IDs saved to: {id_file}")
                    else:
                        print(f"No suitable Dinoflagellate and {target_group} subtrees found in {filename}")
                        
                except Exception as e:
                    print(f"ERROR: Failed to process tree file {filename}: {str(e)}")
                    traceback.print_exc()
                    continue

            summary.write(f"  Processed {processed_trees}/{tree_files_count} tree files\n")
            summary.write(f"  Found {subtrees_found} subtrees\n")
            
            # Extract sequences for this subdirectory
            print(f"\nExtracting sequences for {subdir}...")
            subdir_sequences = extract_dino_sequences(current_dir, default_fasta_file)
            
            # Store sequences for this group if any were found
            if subdir_sequences:
                group_sequences[group_name] = subdir_sequences
                summary.write(f"  Extracted {len(subdir_sequences)} sequences\n")
            else:
                summary.write(f"  No sequences extracted\n")
                
            # Calculate time spent on this directory
            dir_time = time.time() - dir_start_time
            summary.write(f"  Time spent: {dir_time:.2f} seconds\n\n")
            
            # Mark directory as completed
            with open(progress_file, 'a') as pf:
                pf.write(f"Completed processing at {time.ctime()}\n")
                pf.write(f"Time taken: {dir_time:.2f} seconds\n")
                pf.write(f"Processed {processed_trees}/{tree_files_count} tree files\n")
                pf.write(f"Found {subtrees_found} subtrees\n")
                pf.write(f"Extracted {len(subdir_sequences)} sequences\n")

        # Write merged sequences to individual FASTA files
        summary.write("Sequence Merging Results:\n")
        summary.write("----------------------\n")
        
        for group, sequences in group_sequences.items():
            merged_fasta_path = os.path.join(merged_fasta_dir, f"merged_dino_{group}_sequences.fasta")

            with open(merged_fasta_path, 'w') as f:
                for full_id, (header, sequence) in sequences.items():
                    f.write(f">{full_id}\n{sequence}\n")

            print(f"Merged {len(sequences)} Dinoflagellate sequences for {group} into '{merged_fasta_path}'")
            summary.write(f"  {group}: {len(sequences)} sequences\n")

        # Filter merged sequences
        print("\nFiltering merged sequences...")
        summary.write("\nSequence Filtering Results:\n")
        summary.write("-------------------------\n")
        
        for filename in os.listdir(merged_fasta_dir):
            if filename.startswith('merged_dino_') and filename.endswith('_sequences.fasta'):
                input_file = os.path.join(merged_fasta_dir, filename)
                output_file = os.path.join(filtered_fasta_dir, f'filtered_{filename}')

                print(f"\nProcessing {filename}")
                try:
                    unique_combinations, total_sequences, skipped_lines = sort_and_filter_sequences(input_file, output_file)
                    print(f"Processed file: {filename}")
                    print(f"Number of unique branch-OG combinations: {unique_combinations}")
                    print(f"Total number of sequences processed: {total_sequences}")
                    print(f"Number of skipped lines: {skipped_lines}")
                    print(f"Output written to: {output_file}")
                    
                    summary.write(f"  {filename}:\n")
                    summary.write(f"    Unique combinations: {unique_combinations}\n")
                    summary.write(f"    Total sequences: {total_sequences}\n")
                    summary.write(f"    Skipped lines: {skipped_lines}\n")
                except Exception as e:
                    print(f"Error processing {filename}: {str(e)}")
                    summary.write(f"  {filename}: ERROR - {str(e)}\n")

        # Add final summary information
        end_time = time.time()
        total_time = end_time - start_time
        summary.write(f"\nTotal execution time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)\n")
        summary.write(f"Completed at: {time.ctime()}\n")

    print("\nScript execution completed.")
    print(f"Total time: {(time.time() - start_time)/60:.2f} minutes")
    print(f"Summary written to: {summary_file}")
    print(f"Merged sequences: {merged_fasta_dir}")
    print(f"Filtered sequences: {filtered_fasta_dir}")

if __name__ == "__main__":
    main()
