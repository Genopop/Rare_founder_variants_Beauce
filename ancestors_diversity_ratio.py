import geneo as gen
import pandas as pd
import numpy as np
import pickle
import random
import sys
from collections import defaultdict

# === Parameters ===
n_iterations = 1000
n_pro = 7445

# === Input files ===
ped_file = sys.argv[1]
out_file = sys.argv[2]
cache_file = sys.argv[3]

# === Load pedigree ===
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)
print(f"Number of probands: {len(pro)}")

# === Memoized recursive ancestor count ===
def compute_ancestor_counts(ped, individual_id, memo):
    if individual_id in memo:
        return memo[individual_id]

    ancestor_counts = defaultdict(int)
    stack = [(individual_id, False)]
    visited = set()
    temp_results = {}

    while stack:
        current_id, processed = stack.pop()
        if processed:
            try:
                person = ped[int(current_id)]
            except (KeyError, IndexError):
                temp_results[current_id] = defaultdict(int)
                continue

            combined = defaultdict(int)
            for parent in [person.father, person.mother]:
                if parent is not None:
                    pid = parent.ind
                    if pid == 0:
                        continue
                    combined[pid] += 1
                    if pid in temp_results:
                        for anc_id, count in temp_results[pid].items():
                            if anc_id == 0:
                                continue
                            combined[anc_id] += count
            temp_results[current_id] = combined
            continue

        if current_id in memo:
            temp_results[current_id] = memo[current_id]
            continue

        if current_id in visited:
            continue
        visited.add(current_id)

        stack.append((current_id, True))
        try:
            person = ped[int(current_id)]
        except (KeyError, IndexError):
            continue

        for parent in [person.father, person.mother]:
            if parent is not None and parent.ind != 0:
                stack.append((parent.ind, False))

    result = temp_results.get(individual_id, defaultdict(int))
    memo[individual_id] = result
    return result

# === Global memo cache ===
memo = {}

# === Stats: store only cumulative counts to save memory ===
ancestor_stats = defaultdict(lambda: [float('inf'), float('-inf'), 0.0, 0])  # [min, max, sum, count]
iteration_matrix = defaultdict(list)

for i in range(n_iterations):
    print(f"[INFO] Iteration {i+1}/{n_iterations}...")
    sample_pro = random.choices(pro, k=n_pro)

    iteration_counts = defaultdict(int)
    for pro_id in sample_pro:
        anc_counts = compute_ancestor_counts(ped, pro_id, memo)
        for anc_id, count in anc_counts.items():
            if anc_id == 0:
                continue
            iteration_counts[anc_id] += count

    for anc_id, count in iteration_counts.items():
        stats = ancestor_stats[anc_id]
        stats[0] = min(stats[0], count)
        stats[1] = max(stats[1], count)
        stats[2] += count
        stats[3] += 1

        # Update iteration matrix (include zero for missing ancestors)
    all_anc_ids = set(iteration_counts.keys()) | set(iteration_matrix.keys())
    for anc_id in all_anc_ids:
        count = iteration_counts.get(anc_id, 0)
        iteration_matrix[anc_id].append(count)
        
# === Save iteration matrix to CSV ===
df = pd.DataFrame.from_dict(iteration_matrix, orient='index')
df.index.name = 'ancestor_id'
df.columns = [f'iteration_{i+1}' for i in range(n_iterations)]

df.to_csv(out_file)
print(f"[INFO] Iteration matrix saved to {out_file}")

