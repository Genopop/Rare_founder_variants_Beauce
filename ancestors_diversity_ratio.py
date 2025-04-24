## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import psutil
import sys
from sys import argv
import random
from collections import defaultdict

n_iterations=1000
n_pro=7445

def print_memory_usage(stage):
    """Print the current memory usage with a custom message."""
    process = psutil.Process()
    mem_info = process.memory_info()
    print(f"[{stage}] Memory usage: {mem_info.rss / 1024**2:.2f} MB")
  
## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)
print(f"Length of pro: {len(pro)}")

## Identify all unique ancestors
unique_ancestors = gen.ancestor(ped, pro)
print(f"Length of unique_ancestors: {len(unique_ancestors)}")

def get_ancestors_with_counts(ped, individual_id, ancestor_counts=None, visited=None):
    """
    Function to count each time an ancestor appear in the pedigree of all ancestors, accounting for multiple apparition in a same pedigree. 
    Returns : a dictionnary {ancestor_id: number of apparitions}
    """
    if ancestor_counts is None:
        ancestor_counts = defaultdict(int)
    if visited is None:
        visited = set()
    try:
        person = ped[int(individual_id)]
    except (KeyError, IndexError):
        return ancestor_counts  # Indiv not found

    # Avoid going back to the same ancestor
    if individual_id in visited:
        return ancestor_counts
    visited.add(individual_id)

    for parent in [person.father, person.mother]:
        if parent is not None:
            parent_id = parent.ind
            ancestor_counts[parent_id] += 1
            get_ancestors_with_counts(ped, parent_id, ancestor_counts, visited)

    return ancestor_counts

ancestor_stats = defaultdict(list)

for i in range(n_iterations):
    print(f"[INFO] Iteration {i+1}/{n_iterations}...")
    sample_pro = random.choices(pro, k=n_pro)
    ancestor_counts = defaultdict(int)
    for pro_id in sample_pro:
        pro_ancestor_counts = get_ancestors_with_counts(ped, pro_id)
        ## Add count for this iteration
        for anc_id, count in pro_ancestor_counts.items():
            ancestor_counts[anc_id] += count
        ## Add results of the iteration to the dict
    for anc_id, count in ancestor_counts.items():
        ancestor_stats[anc_id].append(count)
    print_memory_usage(f"Iteration {i+1} - After counting ancestors")

# Create DataFrame
ancestor_ids = list(ancestor_stats.keys())
min_counts = [min(ancestor_stats[aid]) for aid in ancestor_ids]
max_counts = [max(ancestor_stats[aid]) for aid in ancestor_ids]
mean_counts = [np.mean(ancestor_stats[aid]) for aid in ancestor_ids]

df_ancestor_stats = pd.DataFrame({
    "AncestorID": ancestor_ids,
    "MinCount": min_counts,
    "MaxCount": max_counts,
    "MeanCount": mean_counts
})

df_ancestor_stats.head()
# Save the file
out_file = sys.argv[2]
df_ancestor_stats.to_csv(out_file, sep='\t', index=False)
