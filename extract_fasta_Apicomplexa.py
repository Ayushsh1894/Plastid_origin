from ete3 import Tree
import os
import re
import sys
from collections import defaultdict

def is_target_leaf(leaf_name, target_group):
    # Extract the taxonomy from folder name (e.g., "Apicomplexa-Archaea_ASGARD" -> "Archaea_ASGARD")
    if isinstance(target_group, list):
        return (leaf_name.startswith("Chromerida") or
                any(leaf_name.startswith(group) for group in target_group))
    else:
        return (leaf_name.startswith("Chromerida") or
                leaf_name.startswith(target_group))

def extract_chromerida_subtrees(tree_file, target_group):
    print(f"Processing file: {tree_file}")
    tree = Tree(tree_file)
    print(f"Tree loaded, number of leaves: {len(tree.get_leaves())}")

    def is_valid_subtree(node):
        children = node.get_children()
        chromerida_child = None
        target_child = None
        other_children = []

        for child in children:
            if all(leaf.name.startswith("Chromerida") for leaf in child.iter_leaves()):
                chromerida_child = child
            # Handle both string and list target groups
            elif isinstance(target_group, list):
                if all(any(leaf.name.startswith(group) for group in target_group) for leaf in child.iter_leaves()):
                    target_child = child
            # If target_group is a string, check if it matches
            elif all(leaf.name.startswith(target_group) for leaf in child.iter_leaves()):
                target_child = child
            else:
                other_children.append(child)

        if chromerida_child and target_child and len(other_children) == 0:
            return True
        elif chromerida_child and target_child and all(is_valid_subtree(child) for child in other_children):
            return True
        return False

    subtrees = []
    for node in tree.traverse(strategy="postorder"):
        if is_valid_subtree(node):
            subtrees.append(node)
            print(f"Found a suitable subtree with {len(node.get_leaves())} leaves")

    print(f"Total suitable subtrees found: {len(subtrees)}")
    return subtrees

def read_fasta(filename):
    sequences = {}
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
    return sequences

def extract_chromerida_sequences(sorted_trees_dir, fasta_file):
    chromerida_sequences = {}
    sequences = read_fasta(fasta_file)

    # Collect all sequence ID files
    for filename in os.listdir(sorted_trees_dir):
        if filename.endswith("_sequence_ids.txt"):
            with open(os.path.join(sorted_trees_dir, filename), 'r') as f:
                for line in f:
                    line = line.strip()
                    og_match = re.search(r'(OG\d+)$', line)
                    if og_match:
                        og_number = og_match.group(1)
                        short_id = line.split('-')[-1].split('_b')[0]  # Extract "Cvel_XXXXX_tX"

                        for header, sequence in sequences.items():
                            if short_id in header:
                                chromerida_sequences[line] = (header, sequence)
                                break

    return chromerida_sequences

def parse_id(id_string):
    # Extract parts from the FASTA header
    match = re.match(r'>Chromerida-(.+)_b(\d+)_OG(\d+)', id_string)
    if match:
        gene_name = match.group(1)
        branch = match.group(2)
        og = match.group(3)
        return gene_name, branch, og
    return None, None, None

def sort_and_filter_sequences(input_file, output_file):
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

def process_fasta_by_organism(input_file, output_dir):
    # Extract the group name from the filename
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    group_name = base_name.split('_')[-2]  # Extract the group name
    
    # Create output files
    os.makedirs(output_dir, exist_ok=True)
    output_cvel = os.path.join(output_dir, f"Api_{group_name}_Cvel.fasta")
    output_vbra = os.path.join(output_dir, f"Api_{group_name}_Vbra.fasta")
    
    with open(input_file, 'r') as infile, open(output_cvel, 'w') as cvel_file, open(output_vbra, 'w') as vbra_file:
        for line in infile:
            if line.startswith('>'):
                # Remove "Chromerida-" from the ID
                line = line.replace("Chromerida-", "")
                # Check if it's a Cvel ID
                if line.startswith('>Cvel'):
                    cvel_file.write(line)
                    current_file = cvel_file
                else:
                    vbra_file.write(line)
                    current_file = vbra_file
            else:
                current_file.write(line)
    
    print(f"Processed {input_file}")
    print(f"Created {output_cvel} and {output_vbra}")
    return output_cvel, output_vbra

def main():
    # Default paths
    base_dir = os.getcwd()  # Current directory (overlaps)
    
    # Check if path is provided as command line argument
    if len(sys.argv) > 1:
        default_fasta_file = sys.argv[1]
    else:
        # Try to read from stdin (for sbatch heredoc)
        try:
            default_fasta_file = input().strip()
        except EOFError:
            default_fasta_file = input("Enter the path to your Chromeride.faa file: ")
            
    print(f"Using FASTA file: {default_fasta_file}")
    
    # Create directories for merged and filtered FASTA files
    merged_fasta_dir = os.path.join(base_dir, "merged_fasta")
    filtered_fasta_dir = os.path.join(base_dir, "filtered_fasta")
    split_fasta_dir = os.path.join(base_dir, "split_fasta_files")
    
    os.makedirs(merged_fasta_dir, exist_ok=True)
    os.makedirs(filtered_fasta_dir, exist_ok=True)
    os.makedirs(split_fasta_dir, exist_ok=True)

    # Create a dictionary to group sequences by their group
    group_sequences = {}

    # Process each directory (Apicomplexa-*)
    for subdir in os.listdir(base_dir):
        if subdir.startswith("Apicomplexa-"):
            current_dir = os.path.join(base_dir, subdir)
            
            # Skip if not a directory
            if not os.path.isdir(current_dir):
                continue
                
            # Get the target group from directory name
            target_group = subdir.split("-")[1]
            
            # Special case for Archaea_ASGARD
            if target_group == "Archaea_ASGARD":
                target_group = ["Archaea", "ASGARD"]
            
            # For organizational purposes
            group_name = target_group.lower() if not isinstance(target_group, list) else "archaea_asgard"
            
            print(f"\nProcessing directory: {subdir}")
            print(f"Target group: {target_group}")

            # Process each tree file in the directory
            tree_files_found = False
            for filename in os.listdir(current_dir):
                if filename.endswith(".treefile"):
                    tree_files_found = True
                    print(f"\nProcessing file: {filename}")
                    tree_file = os.path.join(current_dir, filename)
                    subtrees = extract_chromerida_subtrees(tree_file, target_group)

                    # Extract OG number from filename
                    og_match = re.search(r'(OG\d+)', filename)
                    og_number = og_match.group(1) if og_match else "UnknownOG"

                    if subtrees:
                        print(f"Subtrees for {filename} containing Chromerida and {target_group}:")
                        for i, subtree in enumerate(subtrees, 1):
                            sequence_ids = [leaf.name for leaf in subtree.iter_leaves() if is_target_leaf(leaf.name, target_group)]
                            print(f"Subtree {i}:")
                            print(sequence_ids)
                            print()

                            # Save the subtree
                            outfile = os.path.join(current_dir, f"{og_number}_chromerida_{str(target_group).lower()}_subtree_{i}.nw")
                            subtree.write(outfile=outfile, format=1)
                            print(f"Subtree saved to: {outfile}")

                            # Save sequence IDs to a file
                            id_file = os.path.join(current_dir, f"{og_number}_b{i}_chromerida_{str(target_group).lower()}_sequence_ids.txt")
                            with open(id_file, "w") as f:
                                for seq_id in sequence_ids:
                                    f.write(f"{seq_id}_b{i}_{og_number}\n")
                            print(f"Sequence IDs saved to: {id_file}")
                    else:
                        print(f"No suitable Chromerida and {target_group} subtrees found in {filename}")
            
            if not tree_files_found:
                print(f"No .treefile files found in {current_dir}")
                continue
                
            # Extract sequences for this subdirectory
            subdir_sequences = extract_chromerida_sequences(current_dir, default_fasta_file)
            
            # Store sequences for this group
            group_sequences[group_name] = subdir_sequences

    # Write merged sequences to individual FASTA files
    for group, sequences in group_sequences.items():
        merged_fasta_path = os.path.join(merged_fasta_dir, f"merged_chromerida_{group}_sequences.fasta")

        with open(merged_fasta_path, 'w') as f:
            for full_id, (header, sequence) in sequences.items():
                f.write(f">{full_id}\n{sequence}\n")

        print(f"Merged {len(sequences)} Chromerida sequences for {group} into '{merged_fasta_path}'")

    # Filter merged sequences
    print("\nFiltering merged sequences...")
    for filename in os.listdir(merged_fasta_dir):
        if filename.startswith('merged_chromerida_') and filename.endswith('_sequences.fasta'):
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
            except Exception as e:
                print(f"Error processing {filename}: {e}")
    
    # Process filtered FASTA files to split by organism
    print("\nSplitting filtered sequences by organism...")
    for filename in os.listdir(filtered_fasta_dir):
        if filename.startswith('filtered_merged_chromerida_') and filename.endswith('_sequences.fasta'):
            input_file = os.path.join(filtered_fasta_dir, filename)
            process_fasta_by_organism(input_file, split_fasta_dir)

    print("\nScript execution completed.")
    print(f"Merged sequences: {merged_fasta_dir}")
    print(f"Filtered sequences: {filtered_fasta_dir}")
    print(f"Organism-specific sequences: {split_fasta_dir}")

if __name__ == "__main__":
    main()
