#!/usr/bin/env python3
import os
import sys
import datetime
import shutil
from collections import defaultdict

# Import necessary modules
try:
    from Bio import SeqIO
    from ete3 import Tree
except ImportError:
    print("Required libraries Biopython and ete3 are not installed.")
    print("You can install them using the commands: \n  pip install biopython ete3")
    sys.exit(1)  # Exit if libraries are missing


def is_pure_clade(node, target):
    """Check if a node contains only leaves from a specific group and no other groups."""
    if target == "Archaea_ASGARD":
        target_match = all(leaf.name.startswith("Archaea-") or leaf.name.startswith("ASGARD-") for leaf in node.iter_leaves())
    elif target == "Archaea":
        target_match = all(leaf.name.startswith("Archaea") for leaf in node.iter_leaves())
    elif target == "ASGARD":
        target_match = all(leaf.name.startswith("ASGARD") for leaf in node.iter_leaves()) #seprate check for Archaea and ASGARD
    else:
        target_match = all(leaf.name.startswith(target) for leaf in node.iter_leaves())

    # Also verify the node has at least one leaf (to avoid empty nodes)
    has_leaves = len(list(node.iter_leaves())) > 0

    return target_match and has_leaves


def has_both_groups(node, group1, group2):
    """Check if a node contains sequences from both specified groups."""
    has_group1 = any(leaf.name.startswith(group1) for leaf in node.iter_leaves())
    if group2 == "Archaea_ASGARD":
        has_group2 = any(leaf.name.startswith("Archaea-") or leaf.name.startswith("ASGARD-") for leaf in node.iter_leaves())
    elif group2 == "Archaea":
         has_group2 = any(leaf.name.startswith("Archaea") for leaf in node.iter_leaves())
    elif group2 == "ASGARD":
         has_group2 = any(leaf.name.startswith("ASGARD") for leaf in node.iter_leaves()) #seprate check for Archaea and ASGARD
    else:
        has_group2 = any(leaf.name.startswith(group2) for leaf in node.iter_leaves())
    return has_group1 and has_group2


def find_mixed_pure_sister_clades(tree, main_target_group):
    """Find where a (Dinoflagellates + main_target_group) is sister to a pure main_target_group."""
    results = []

    for node in tree.traverse("postorder"):
        if len(node.children) != 2:  # Only consider nodes with exactly two children
            continue

        child1, child2 = node.children

        # Verify mixed node has ONLY Dinoflagellates and target group (no other taxa)
        def is_valid_mixed_node(node, group1, group2):
            # Check that every leaf in the node is either from group1 or group2 (no other groups)
            for leaf in node.iter_leaves():
                if group2 == "Archaea_ASGARD":
                    is_group2 = leaf.name.startswith("Archaea-") or leaf.name.startswith("ASGARD-")
                elif group2 == "Archaea":
                    is_group2 = leaf.name.startswith("Archaea")
                elif group2 == "ASGARD":
                    is_group2 = leaf.name.startswith("ASGARD") #seprate check for Archaea and ASGARD
                else:
                    is_group2 = leaf.name.startswith(group2)

                if not (leaf.name.startswith(group1) or is_group2):
                    return False

            # Also verify it has both groups as required
            return has_both_groups(node, group1, group2)

        # Case 1: child1 is STRICTLY (Dinoflagellates + target), child2 is pure target
        if (is_valid_mixed_node(child1, "Dinoflagellates", main_target_group) and
            is_pure_clade(child2, main_target_group)):
            results.append((node, child1, child2))

        # Case 2: child2 is STRICTLY (Dinoflagellates + target), child1 is pure target
        elif (is_valid_mixed_node(child2, "Dinoflagellates", main_target_group) and
              is_pure_clade(child1, main_target_group)):
            results.append((node, child2, child1))

    return results


def find_dinoflagellates_targetgroup_sister_to_pure_dinoflagellates(tree, target_group):
    """Find where a (Dinoflagellates + specified target group) is sister to a pure Dinoflagellates clade."""
    results = []

    for node in tree.traverse("postorder"):
        if len(node.children) != 2:  # Only consider nodes with exactly two children
            continue

        child1, child2 = node.children

        # Verify mixed node has ONLY Dinoflagellates and the target group (no other taxa)
        def is_valid_mixed_node(node, target_group):
            # Check that every leaf in the node is either Dinoflagellates or the target group
            for leaf in node.iter_leaves():
                is_target = False
                if target_group == "Archaea_ASGARD":
                    is_target = leaf.name.startswith("Archaea-") or leaf.name.startswith("ASGARD-")
                elif target_group == "Archaea":
                    is_target = leaf.name.startswith("Archaea")
                elif target_group == "ASGARD":
                    is_target = leaf.name.startswith("ASGARD")
                else:
                    is_target = leaf.name.startswith(target_group)
                
                if not (leaf.name.startswith("Dinoflagellates") or is_target):
                    return False

            # Also verify it has both groups as required
            return has_both_groups(node, "Dinoflagellates", target_group)

        # Case 1: child1 is STRICTLY (Dinoflagellates + target), child2 is pure Dinoflagellates
        if (is_valid_mixed_node(child1, target_group) and 
            is_pure_clade(child2, "Dinoflagellates")):
            results.append((node, child1, child2))

        # Case 2: child2 is STRICTLY (Dinoflagellates + target), child1 is pure Dinoflagellates
        elif (is_valid_mixed_node(child2, target_group) and 
              is_pure_clade(child1, "Dinoflagellates")):
            results.append((node, child2, child1))

    return results


def extract_dinoflagellates_ids(node):
    """Extract Dinoflagellates leaf names from a node."""
    return [leaf.name for leaf in node.iter_leaves() if leaf.name.startswith("Dinoflagellates")]


def get_target_group_from_folder(folder_name):
    """Extract the target group from the folder name (after the hyphen)."""
    if "-" not in folder_name:
        return None
    
    return folder_name.split("-", 1)[1]


def process_single_folder(input_dir, output_base_dir, target_groups, dinoflagellates_fasta_file):
    """Process a single input folder and produce separate results."""
    
    folder_name = os.path.basename(input_dir)
    print(f"\n--- Processing folder: {folder_name} ---")
    
    # Extract the target group for this folder (the part after "Dinoflagellates-")
    folder_target_group = get_target_group_from_folder(folder_name)
    if not folder_target_group:
        print(f"Warning: Could not determine target group from folder name {folder_name}")
        folder_target_group = "Unknown"
    else:
        print(f"Target group extracted from folder name: {folder_target_group}")
    
    if not os.path.isdir(input_dir):
        print(f"Skipping - Input directory not found: {input_dir}")
        return {}, {}, 0
    
    # Create folder-specific output directory
    folder_output_dir = os.path.join(output_base_dir, folder_name)
    os.makedirs(folder_output_dir, exist_ok=True)

    # Initialize tracking variables
    all_results = {}
    all_dinoflagellates_ids = {}
    
    # Special variables for Dinoflagellates-TargetGroup relationship
    dinoflagellates_trees_count = 0
    dinoflagellates_sequences = set()
    
    # Special variables for merged Archaea+ASGARD
    merged_archaea_asgard_trees = set()  # Use a set to avoid duplicates
    merged_archaea_asgard_sequences = set()
    
    # Files processed count
    folder_files_count = 0
    
    # Create special output directory for Dinoflagellates
    dinoflagellates_output_dir = os.path.join(folder_output_dir, "Valid_Dinoflagellates_Sister_Trees")
    os.makedirs(dinoflagellates_output_dir, exist_ok=True)
    
    # Create merged output directory for Archaea+ASGARD if applicable
    if folder_target_group == "Archaea_ASGARD":
        merged_archaea_asgard_dir = os.path.join(folder_output_dir, "Valid_Archaea_ASGARD_Merged_Sister_Trees")
        os.makedirs(merged_archaea_asgard_dir, exist_ok=True)
    
    # Process all tree files in the directory
    for filename in os.listdir(input_dir):
        if not filename.endswith(".treefile"):
            continue  # Skip non-tree files
            
        folder_files_count += 1
        tree_file = os.path.join(input_dir, filename)

        try:
            print(f"Processing tree: {filename}")
            tree = Tree(tree_file, format=1)
            
            # Check for the specific (Dinoflagellates + folder_target_group) sister to pure Dinoflagellates case
            dinoflagellates_nodes = find_dinoflagellates_targetgroup_sister_to_pure_dinoflagellates(tree, folder_target_group)
            if dinoflagellates_nodes:
                # Copy valid tree to the Dinoflagellates output directory
                shutil.copy2(tree_file, os.path.join(dinoflagellates_output_dir, filename))
                dinoflagellates_trees_count += 1
                
                # Extract Dinoflagellates IDs
                for mixed_child, _, _ in dinoflagellates_nodes:
                    dinoflagellates_sequences.update(extract_dinoflagellates_ids(mixed_child))
                
                print(f"  Found valid (Dinoflagellates + {folder_target_group}) sister to (pure Dinoflagellates) relationship.")
            
            # Process all target groups to find (Dinoflagellates + target) sister to pure target
            for main_target_group in target_groups:
                # Skip Dinoflagellates as we're handling it separately
                # Skip Archaea_ASGARD as we're not using it anymore
                if main_target_group in ["Dinoflagellates", "Archaea_ASGARD"]:
                    continue
                    
                # Initialize for this target group if not already done
                if main_target_group not in all_results:
                    all_results[main_target_group] = 0
                    all_dinoflagellates_ids[main_target_group] = set()

                # Create an output directory for valid main target group trees
                group_output_dir = os.path.join(folder_output_dir, f"Valid_{main_target_group}_Sister_Trees")
                os.makedirs(group_output_dir, exist_ok=True)

                # Find patterns (Dinoflagellates + main_target_group) sister to (pure main_target_group)
                valid_nodes = find_mixed_pure_sister_clades(tree, main_target_group)

                if valid_nodes:
                    # Copy valid tree to output directory
                    shutil.copy2(tree_file, os.path.join(group_output_dir, filename))
                    all_results[main_target_group] += 1

                    # Extract Dinoflagellates IDs and add to the set for this group
                    for _, mixed_child, _ in valid_nodes:
                        all_dinoflagellates_ids[main_target_group].update(extract_dinoflagellates_ids(mixed_child))

                    print(f"  Found valid (Dinoflagellates + {main_target_group}) sister to (pure {main_target_group}) relationship.")
                    
                    # Special handling for Archaea and ASGARD - merge them
                    if main_target_group in ["Archaea", "ASGARD"] and folder_target_group == "Archaea_ASGARD":
                        merged_archaea_asgard_trees.add(tree_file)
                        extracted_ids = extract_dinoflagellates_ids(mixed_child)
                        merged_archaea_asgard_sequences.update(extracted_ids)

        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Add Dinoflagellates results to the main results
    all_results["Dinoflagellates"] = dinoflagellates_trees_count
    all_dinoflagellates_ids["Dinoflagellates"] = dinoflagellates_sequences
    
    # Add merged Archaea+ASGARD results if applicable
    if folder_target_group == "Archaea_ASGARD":
        all_results["Archaea_ASGARD_Merged"] = len(merged_archaea_asgard_trees)
        all_dinoflagellates_ids["Archaea_ASGARD_Merged"] = merged_archaea_asgard_sequences
        
        # Copy trees to the merged directory
        for tree_file in merged_archaea_asgard_trees:
            shutil.copy2(tree_file, os.path.join(merged_archaea_asgard_dir, os.path.basename(tree_file)))

    # Extract Dinoflagellates sequences for each group
    for main_target_group in list(all_results.keys()):
        if not all_dinoflagellates_ids[main_target_group]:
            print(f"No Dinoflagellates sequences found for {main_target_group}. Skipping sequence extraction.")
            continue
            
        print(f"Extracting Dinoflagellates sequences for {main_target_group}...")
        extracted_sequences = []
        with open(dinoflagellates_fasta_file, "r") as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                if record.id in all_dinoflagellates_ids[main_target_group]:
                    extracted_sequences.append(record)

        output_fasta = os.path.join(folder_output_dir, f"Extracted_Dinoflagellates_{main_target_group}_sequences.faa")
        SeqIO.write(extracted_sequences, output_fasta, "fasta")
        print(f"Extracted {len(extracted_sequences)} Dinoflagellates sequences to {output_fasta}.")
    
    # Create summary report for this folder
    summary_report(all_results, all_dinoflagellates_ids, folder_files_count, input_dir, folder_output_dir, folder_target_group)
    
    print(f"Processed {folder_files_count} tree files in {folder_name}")
    
    return all_results, all_dinoflagellates_ids, folder_files_count


def summary_report(results, dinoflagellates_ids, total_files, input_dir, output_dir, folder_target_group):
    """Generate a summary report for sister relationships."""
    summary_file = os.path.join(output_dir, "Sister_Relationships_Summary.txt")
    with open(summary_file, "w") as f:
        f.write("== Sister Relationships Analysis Summary ==\n")
        f.write("Generated on: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
        
        # Write input folder processed
        f.write("== Input Folder ==\n")
        f.write(f"- {os.path.basename(input_dir)}\n\n")
        
        f.write(f"Total tree files processed: {total_files}\n\n")
        
        # Write summary for each target group
        f.write("== Results by Target Group ==\n")
        
        # Sort groups by count for better readability
        sorted_groups = sorted(results.items(), key=lambda x: x[1], reverse=True)
        
        for group, count in sorted_groups:
            # Skip groups with no hits
            if count == 0:
                continue
                
            f.write(f"{group}:\n")
            f.write(f"  Trees found: {count}\n")
            
            if group == "Dinoflagellates":
                f.write(f"  Unique Dinoflagellates sequences: {len(dinoflagellates_ids[group])}\n")
                f.write(f"  Note: These are cases where (Dinoflagellates + {folder_target_group}) clades are sister to pure Dinoflagellates clades.\n")
            elif group == "Archaea_ASGARD_Merged":
                f.write(f"  Unique Dinoflagellates sequences: {len(dinoflagellates_ids[group])}\n")
                f.write(f"  Note: This represents the merged results from both Archaea and ASGARD groups.\n")
            else:
                f.write(f"  Unique Dinoflagellates sequences: {len(dinoflagellates_ids[group])}\n")
            
            f.write("\n")
        
        # Add overall statistics
        total_trees = sum(results.values())
        # Adjust for double-counting of Archaea and ASGARD in the merged group
        if "Archaea" in results and "ASGARD" in results and "Archaea_ASGARD_Merged" in results:
            total_trees -= results["Archaea_ASGARD_Merged"]
            f.write("Note: Total adjusted to avoid double-counting trees in the Archaea_ASGARD_Merged category.\n")
            
        f.write("== Overall Statistics ==\n")
        f.write(f"Total trees with valid relationships: {total_trees}\n")
        
        # Calculate total unique Dinoflagellates sequences across all groups
        all_dinoflagellates = set()
        for group, ids in dinoflagellates_ids.items():
            # Skip Archaea_ASGARD_Merged to avoid double-counting
            if group != "Archaea_ASGARD_Merged":
                all_dinoflagellates.update(ids)
        
        f.write(f"Total unique Dinoflagellates sequences: {len(all_dinoflagellates)}\n")
    
    print(f"Summary report created: {summary_file}")


def main():
    """Main function to run the analysis."""
    # Define target groups to analyze
    target_groups = [
        "Dinoflagellates",
        "Rhizaria", 
        "Stramenopiles", 
        "Glaucocystophyceae", 
        "Cryptophyta", 
        "Haptophyta", 
        "Rhodophyta", 
        "Discoba", 
        "Rhodelphis", 
        "Opisthokonta", 
        "Metamonada", 
        "Archaea", 
        "ASGARD", 
        "Archaea_ASGARD",
        "Evosea", 
        "Proteobacteria", 
        "Viridiplantae", 
        "Cyanobacteria",
        "Chromerida"
    ]
    
    # Define paths
    base_dir = os.getcwd()
    
    # Define input folders (all Dinoflagellates-* folders)
    input_folders = [
        os.path.join(base_dir, "Dinoflagellates-Archaea_ASGARD"),
        os.path.join(base_dir, "Dinoflagellates-Cryptophyta"),
        os.path.join(base_dir, "Dinoflagellates-Cyanobacteria"),
        os.path.join(base_dir, "Dinoflagellates-Discoba"),
        os.path.join(base_dir, "Dinoflagellates-Evosea"),
        os.path.join(base_dir, "Dinoflagellates-Glaucocystophyceae"),
        os.path.join(base_dir, "Dinoflagellates-Haptophyta"),
        os.path.join(base_dir, "Dinoflagellates-Metamonada"),
        os.path.join(base_dir, "Dinoflagellates-Opisthokonta"),
        os.path.join(base_dir, "Dinoflagellates-Proteobacteria"),
        os.path.join(base_dir, "Dinoflagellates-Rhizaria"),
        os.path.join(base_dir, "Dinoflagellates-Rhodelphis"),
        os.path.join(base_dir, "Dinoflagellates-Rhodophyta"),
        os.path.join(base_dir, "Dinoflagellates-Stramenopiles"),
        os.path.join(base_dir, "Dinoflagellates-Viridiplantae")
    ]
    
    output_base_dir = os.path.join(base_dir, "Separate_Dinoflagellates_Analysis_Results")
    dinoflagellates_fasta_file = "/home/users/ayush/Ayush/New/DBs/Dinoflagellate.faa"  # Path to the dinoflagellates FASTA file
    
    # Create main output directory
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Validate FASTA file path
    if not os.path.isfile(dinoflagellates_fasta_file):
        print(f"Dinoflagellates FASTA file not found: {dinoflagellates_fasta_file}")
        sys.exit(1)
    
    # Process each folder separately
    total_files_processed = 0
    
    for input_dir in input_folders:
        if os.path.isdir(input_dir):
            _, _, folder_files_count = process_single_folder(input_dir, output_base_dir, target_groups, dinoflagellates_fasta_file)
            total_files_processed += folder_files_count
    
    print(f"\nAnalysis completed successfully!")
    print(f"Total files processed: {total_files_processed}")
    print(f"Results saved to: {output_base_dir}")


# Run the main function when the script is executed
if __name__ == "__main__":
    main()
