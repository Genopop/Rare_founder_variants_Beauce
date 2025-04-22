## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

print("2nd step: Filter for completeness")
## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
pro = gen.pro(ped)

## Get completeness at desired generation
gen_threshold = int(sys.argv[2])
comp = gen.completeness(ped, pro=pro, genNo=gen_threshold, type="IND")
print(f"Length of completeness: ", len(comp))
df_comp = pd.DataFrame({
    'proID': pro,
    'compl': comp
})
print(df_comp.head())

## Keep complete individuals only
comp_threshold = int(sys.argv[3])
df_comp_filtered = df_comp[df_comp['compl'] >= comp_threshold]
comp_pro = df_comp_filtered['proID']
print(f"Length of probands with sufficient completeness: {len(comp_pro)}")

## Branch the genealogy
new_ped = gen.branching(ped, pro=comp_pro)
new_asc = gen.genout(new_ped)

## Save the new asc
print(new_asc.head())
inb_file = sys.argv[4]
new_asc.to_csv(
        inb_file,
        sep='\t',
        index=False
)

## Compute mean completeness (for the plot)
mean_file = sys.argv[5]
mean_compl = gen.completeness(new_ped, type = "MEAN")
mean_compl_df = pd.DataFrame(mean_compl)
mean_compl_df.to_csv(
        mean_file,
        sep='\t',
        index=False
)
