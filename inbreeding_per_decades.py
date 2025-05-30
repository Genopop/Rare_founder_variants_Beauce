## Import packages
import geneakit as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)

## Get all IDs and convert to list
IDs = list(ped.keys())
print(f"Length of IDs: {len(IDs)}")
anc = gen.ancestor(ped, pro)
anc_df = pd.DataFrame({'IndID': anc}).astype({'IndID': str})
print("Length of ancestors: ", len(anc))

## Add marriage date infos
df = pd.read_csv("dates_locations.txt",
	delimiter=' ', dtype={'IndID': str})  
anc_info = anc_df.merge(df, on='IndID')
## Add decades and round to the nearest upper bound
anc_info['decades'] = np.ceil(anc_info['DateMariage'] / 10) * 10
decades = anc_info['decades'].unique()
## Filter out individuals with decades == 0 or decades > 1960
anc_info = anc_info[(anc_info['decades'] != 0) & (anc_info['decades'] <= 1960)]
## List ancestors per decades
grouped = anc_info.groupby('decades')
## Create a list of IndID lists, one per decade
decade_data = [(decade, group['IndID'].tolist()) for decade, group in grouped]

## Compute inbreeding and confidence interval for each decade
meanf_per_decade = []
CI_per_decade = []
print("Inbreeding computation...")
for decade, IDs in decade_data:
    print(f"Processing Decade {decade} with {len(IDs)} individuals...")
    f = gen.f(ped, pro=IDs)
    f_values = f["F"].astype(float).tolist()
    max_f = max(f_values)
    print(f"Maximum value: {max_f}")
    meanf = sum(f_values) / len(f_values)
    print(f"Mean value: {meanf}")
    meanf_per_decade.append({"decade": decade, "mean": meanf})
    if max_f == 0.0:
        CIf = np.nan  # Assign NaN instead of computing if max(f_values) is 0
    else:
        CIf = gen.fCI(f)
    CI_per_decade.append({"decade": decade, "CI": CIf})

## Save as CSV
meanf_df = pd.DataFrame(meanf_per_decade)
print(meanf_df.head())
inb_file = sys.argv[2]
meanf_df.to_csv(
        inb_file,
        sep=',',
        index=False
)

data_to_save = pd.DataFrame(CI_per_decade)
CI_file = sys.argv[3]
data_to_save.to_csv(
	CI_file,
	sep=",",
	index=False
)

