## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

print("1st step: Subset the pedigree file into regional groups")
## Prepare the whole file to subset
whole_ped = gen.genealogy("/path/to/all_pedigrees.asc")
all_pro = gen.pro(whole_ped)

## Region infos
regions = pd.read_csv(
    "/patho/to/all_pedigrees_dates_locations.txt", 
    sep=' ',
    dtype={'RegionIDMariage': int, 'IndID': int}
)

## Identify probands from the same region 
regional_code=sys.argv[1]
regional_ids = regions[regions['RegionIDMariage'] == int(regional_code)]
regional_probands = set(all_pro).intersection(set(regional_ids['IndID']))
regional_probands = list(regional_probands)
print(f"Length of probands married in the regional subset: {len(regional_probands)}")

## Save the new pedigree information
regional_ped = gen.branching(whole_ped, pro=regional_probands)
regional_asc = gen.genout(regional_ped)
out_file = sys.argv[2]
regional_asc.to_csv(
        out_file,
        sep='\t',
        index=False
)
