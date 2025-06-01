import pandas as pd
import os
from collections import defaultdict

# === Input: List of IBD file paths ===
file_list_path = "ibd_files_list.txt"
# === Output files ===
outfile_merged = os.path.join(path_out, "all_chr.merged.ibd.txt")
outfile_proportion = os.path.join(path_out, "all_chr.merged.ibd_with_prop_seg_shared.txt")
# === Initialization ===
header_written = False
ibd_totals = defaultdict(float)
# === Read file list ===
with open(file_list_path, "r") as f:
    file_names = [line.strip() for line in f if line.strip()]
# === Process each IBD file ===
for file in file_names:
    print(f"Processing file: {file}")
    try:
        # Read IBD data with expected columns
        df = pd.read_csv(
            file,
            sep="\t",
            header=0,
            names=['ind1', 'chr1', 'ind2', 'chr2', 'chr', 'start', 'end', 'LOD', 'length'],
            dtype=str
        )
        # Replace decimal separator and convert to numeric
        df['length'] = df['length'].str.replace(",", ".", regex=False)
        df['length'] = pd.to_numeric(df['length'], errors='coerce')
        # Filter segments with length >= 2 cM
        df = df[df['length'] >= 2]
        # Debug preview
        print(df.head())
        # Append to merged IBD file
        df.to_csv(outfile_merged, sep="\t", index=False, mode='a', header=not header_written)
        header_written = True
        # Accumulate total IBD per pair
        for _, row in df.iterrows():
            key = (row['ind1'], row['ind2'])
            ibd_totals[key] += row['length']
    except Exception as e:
        print(f"Error reading file: {file}: {e}")

# === Compute summary statistics ===
print("Generating summary file")
length_genome = max(ibd_totals.values())
print(length_genome)
# Build DataFrame of total shared segment lengths
total_IBD_df = pd.DataFrame(
    [(k[0], k[1], v) for k, v in ibd_totals.items()],
    columns=['ind1', 'ind2', 'total_length']
)
# Calculate percentage of genome shared
total_IBD_df["%seg_shared"] = (total_IBD_df["total_length"] / length_genome) * 100
# Save result
total_IBD_df.to_csv(outfile_proportion, sep="\t", index=False)

print("Process completed successfully.")
