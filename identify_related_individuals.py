# Import necessary packages
import geneo as gen
import pandas as pd
import numpy as np
import pickle
import networkx as nx
import csv

# Load pedigree data
pedigree_path = "/lustre03/project/6033529/quebec_10x/data/genealogy/gen.Janv2022.aveCodesGenetiques.csv"
pedigree_df = pd.read_csv(pedigree_path, sep=";")
pedigree = gen.genealogy(pedigree_df)
# Load list of genotyped individuals
ids_path = "/lustre03/project/6033529/schizo/mylgag/Beauce_founder_effect/enriched_variants/data/schizo.cag.hg38.commonsnps.mind0.05.hwe10-6_LDpruned2_maf0.05.fam"
ids_df = pd.read_csv(ids_path, sep="\s+", header=None)
genotyped_ids = ids_df.iloc[:, 1]
# Filter pedigree for genotyped individuals
genotyped_df = pedigree_df[pedigree_df['ind'].isin(genotyped_ids)]
# Create subset of the pedigree
subset_pedigree = gen.branching(pedigree, pro=genotyped_df['ind'])

# Compute kinship matrix and save it
kinship_matrix = gen.phi(subset_pedigree)
kinship_output_path = "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/genealogies/results/schizo_kinship_common_inds.pkl"
with open(kinship_output_path, "wb") as f:
    pickle.dump(kinship_matrix, f)

# Build kinship graph (threshold: 0.0625)
kinship_graph = nx.Graph()
kinship_graph.add_nodes_from(kinship_matrix.index)
for i in range(len(kinship_matrix)):
    for j in range(i + 1, len(kinship_matrix)):
        if kinship_matrix.iloc[i, j] >= 0.0625:
            kinship_graph.add_edge(kinship_matrix.index[i], kinship_matrix.columns[j])
# Extract kinship-based groups
kinship_groups = list(nx.connected_components(kinship_graph))

# Load IBD results
ibd_path = "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/results_clinvar/all_chr.merged.ibd_with_prop_seg_shared.txt"
ibd_df = pd.read_csv(ibd_path, sep="\t")
ibd_df['%seg_shared'] = pd.to_numeric(ibd_df['%seg_shared'], errors='coerce')
third_degree_ibd = ibd_df[ibd_df['%seg_shared'] >= 12.5]

# Build IBD graph
ibd_graph = nx.Graph()
for _, row in third_degree_ibd.iterrows():
    ibd_graph.add_edge(row['ind1'], row['ind2'])
# Extract IBD-based groups
ibd_groups_raw = list(nx.connected_components(ibd_graph))

# Print IBD-based groups
print("Groups of third-degree related individuals:")
for group in ibd_groups_raw:
    print(group)

# Save IBD groups to file
ibd_groups_output_path = "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/IBD_related_groups.txt"
with open(ibd_groups_output_path, "w") as f:
    for group in ibd_groups_raw:
        f.write(",".join(group) + "\n")

# Load saved IBD groups and normalize IDs (remove suffixes)
ibd_groups = []
with open(ibd_groups_output_path, "r") as f:
    for line in f:
        members = line.strip().split(",")
        group = sorted({m.split("_")[0] for m in members}, key=int)
        ibd_groups.append(group)

# Convert groups to sets for comparison
kinship_group_sets = [set(map(str, group)) for group in kinship_groups]
ibd_group_sets = [set(map(str, group)) for group in ibd_groups]

# Identify overlapping individuals between kinship and IBD groups
intersections = []
for kin_group in kinship_group_sets:
    for ibd_group in ibd_group_sets:
        common_inds = kin_group & ibd_group
        if len(common_inds) >= 2:
            intersections.append(common_inds)

# Evaluate concordance between methods
all_ibd_inds = set().union(*ibd_group_sets)
for idx, kin_group in enumerate(kinship_group_sets, start=1):
    in_ibd = kin_group & all_ibd_inds
    not_in_ibd = kin_group - all_ibd_inds
    print(f"Groupe {idx}:")
    print(f"  - Detected by refinedIBD : {', '.join(sorted(in_ibd, key=int)) if in_ibd else 'Aucun'}")
    print(f"  - Not detected by refinedIBD : {', '.join(sorted(not_in_ibd, key=int)) if not_in_ibd else 'Aucun'}\n")

# Merge kinship and IBD groups, ensuring uniqueness
all_groups_combined = kinship_group_sets + ibd_group_sets
unique_groups = []
seen = set()

for group in all_groups_combined:
    group_frozen = frozenset(group)
    if group_frozen not in seen:
        seen.add(group_frozen)
        unique_groups.append(group)

# Display final groups
for idx, group in enumerate(unique_groups, start=1):
    print(f"Group {idx}: {', '.join(sorted(group, key=int))}")

# Save merged related groups
group_output_path = "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/related_groups_refinedIBD0.125_kinship0.0625.csv"
with open(group_output_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows([sorted(group, key=int) for group in unique_groups])

# Select one individual per group (maximize unrelatedness)
selected_individuals = set()
final_selection = []

for group in unique_groups:
    for individual in sorted(group):
        if individual not in selected_individuals:
            final_selection.append(individual)
            selected_individuals.add(individual)
            break  # Only select one per group

print("Selected individuals (maximizing unique picks):")
print(final_selection)

# Identify individuals to exclude from projection
all_related_inds = set().union(*unique_groups)
excluded_inds = all_related_inds - set(final_selection)
excluded_df = pd.DataFrame({'excluded_individual': sorted(excluded_inds)})
excluded_split_df = excluded_df['excluded_individual'].str.split('_', expand=True)

# Save individuals to exclude
exclusion_output_path = "/home/mylgag/projects/rrg-girardsi/schizo/mylgag/Beauce_founder_effect/enriched_variants/results/IBD_related_groups0.125_and_pedigree_groups0.0625_inds_to_exclude.txt"
excluded_split_df.to_csv(exclusion_output_path, sep='\t', index=False, header=False)
