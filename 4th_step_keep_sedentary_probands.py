## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

print("4th step: Only keep sedentary probands")
## "Sedentary probands" refers to individuals married in the same region as their 2 parents and 4 grand parents ##
## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)

## Get the parents' IDs
father_ids = [ped[individual_id].father for individual_id in pro]
mother_ids = [ped[individual_id].mother for individual_id in pro]
parents_ids = father_ids + mother_ids
print(f"Length of probands' parents: {len(parents_ids)}")


## Ensure that parents were in the same region
regions = pd.read_csv(
    "/path/to/all_pedigrees_dates_locations.txt", 
    sep=' ',
    dtype={'RegionIDMariage': int, 'IndID': int}
)

regional_code=sys.argv[2]
regional_ids = regions[regions['RegionIDMariage'] == int(regional_code)]
# Extract only father and mother IDs (accessing the .ind attribute)
parents_set = {entry.father.ind for entry in parents_ids}.union({entry.mother.ind for entry in parents_ids})
# Find intersection
regional_parents = list(parents_set.intersection(set(regional_ids['IndID'])))
print(f"Length of probands' parents married in the region: {len(regional_parents)}")

## Get the grandparents' IDs
gf_ids = [ped[individual_id].father.ind for individual_id in regional_parents if ped[individual_id].father is not None]
gm_ids = [ped[individual_id].mother.ind for individual_id in regional_parents if ped[individual_id].mother is not None]
gp_ids = gf_ids + gm_ids
print(f"Length of probands' parents married in the region grandparents: {len(gp_ids)}")

## Check if grandparents were in the same region
gp_regional = set(gp_ids).intersection(set(regional_ids['IndID']))
gp_regional = list(gp_regional)  # Convert to list
print(f"Length of probands' grandparents married in the region: {len(gp_regional)}")

valid_probands = []
## Iterate through each proband
for individual_id in pro:
	# Get the parents of the proband
	father = ped[individual_id].father
	mother = ped[individual_id].mother

	# Get the grandparents of the proband
	gf_father = ped[father].father if father in ped else None
	gm_father = ped[father].mother if father in ped else None
	gf_mother = ped[mother].father if mother in ped else None
	gm_mother = ped[mother].mother if mother in ped else None

	# Collect all grandparents' IDs
	grandparents = [gp for gp in [gf_father, gm_father, gf_mother, gm_mother] if gp is not None]
	# Check if all four grandparents are in the gp_regional list
	if all(gp in gp_regional for gp in grandparents):
		valid_probands.append(individual_id)  # Add proband if condition is met

print(f"Number of probands with all 4 grandparents married in region: {len(valid_probands)}")
valid_probands = list(map(int, valid_probands))

## Save the new pedigree information
sed_ped = gen.branching(ped, pro=valid_probands)
sed_asc = gen.genout(sed_ped)
out_file = sys.argv[3]
sed_asc.to_csv(
        out_file,
        sep='\t',
        index=False
)
