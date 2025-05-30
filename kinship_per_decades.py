## Import packages
import geneakit as gen
import pandas as pd
import numpy as np
import pickle
import sys
from sys import argv

## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)
anc = gen.ancestor(ped, pro)
anc_df = pd.DataFrame({'IndID': anc}).astype({'IndID': str}) 
print("Length of ancestors: ", len(anc))

## Add marriage date infos
df = pd.read_csv("dates_locations.txt", 
                 delimiter=' ', dtype={'IndID': str})
anc_info = anc_df.merge(df, on='IndID')
## Add decades and round to the nearest
anc_info['decades'] = np.ceil(anc_info['DateMariage'] / 10) * 10
decades = anc_info['decades'].unique()
## Filter out individuals with decades == 0 or decades > 1960
anc_info = anc_info[(anc_info['decades'] != 0) & (anc_info['decades'] <= 1960)]
## List ancestors per decades
grouped_df = anc_info.groupby('decades')
## Create a list of IndID lists, one per decade
list_of_ind = [group['IndID'].tolist() for _, group in grouped_df]

## Remove the siblings from each slice
final_list = []
print("Siblings removal...")
for inds in list_of_ind:
	tokeep = set()
	print("Length of ancestors: ", len(inds))

	inds = inds.copy()
	while inds:
		ID = inds.pop()
		sibs = gen.sibship(ped, [ID])
		if not tokeep.intersection(sibs):  # Check if no siblings are already in the set
			tokeep.add(ID)

	print("Length of ancestors to keep: ", len(tokeep))
	# Convert to a vector and store the result in final_list
	final_list.append(sorted(tokeep))

## Remove the children from each slice
no_rel_list = []
print("Children removal...")
for inds in final_list:
	tokeep = set()
	inds = inds.copy()
	while inds:
		ID = inds.pop()
		kids = gen.children(ped, [int(ID)])
		if not tokeep.intersection(kids): # Check if no children are already in the set
			tokeep.add(int(ID))
	print("Length of final ancestors to keep: ", len(tokeep))
	# Convert to a vector and store the result in final_list
	no_rel_list.append(sorted(tokeep))

## Prepare a list of unrelated individuals per decade
decades = list(grouped_df.groups.keys())
no_rel_list = [group['IndID'].tolist() for _, group in grouped_df if len(group) > 0]
print(f"Length of decades: {len(decades)}")
print(f"Length of no_rel_list: {len(no_rel_list)}")

decades = [int(d) for d in decades]
output_df = pd.DataFrame({
    'decade': decades,
    'IndID_list': [", ".join(map(str, inds)) for inds in no_rel_list]  # Convert lists to comma-separated strings
})
out_file = sys.argv[2]
output_df.to_csv(f"{out_file}.CSV", index=False)

## Get mean kinship per decade
kin_dir="/path/where/to/save/kinship/matrices/"
mean_phi_per_decade = []
# CI_per_decade = []
print("Mean kinship computation...")
for i, IDs in enumerate(no_rel_list):
	decade = df['decade'][i]
	print(f"Processing Decade {decade} with {len(IDs)} individuals...")
	phi = gen.phi(ped, pro=IDs)
	mean_phi = gen.phiMean(phi)
    
	print(f"Mean kinship for Decade {df['decade'][i]}: {mean_phi}")
	mean_phi_per_decade.append(mean_phi)

	pickle_filename = f"{kin_dir}phi_mtl_decade_{decade}.pkl"
	with open(pickle_filename, "wb") as f:
		pickle.dump(phi, f)
	del phi # Free up memory

## Save mean kinship
data_to_save = pd.DataFrame(mean_phi_per_decade)
kin_file = sys.argv[3]
data_to_save.to_csv(kin_file)
