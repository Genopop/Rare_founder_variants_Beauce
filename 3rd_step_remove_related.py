## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

print("3rd step: Remove probands related to the 1st degree")
## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)

## Identify the non-related probands
tokeep = set()
stack = gen.pro(ped)

while stack:
	ID = stack.pop()
	sibs = gen.sibship(ped, [ID])

	if not tokeep.intersection(sibs):
		tokeep.add(ID)

## Convert set to list
tokeep = list(tokeep)
print(f"Length of unrelated probands to keep: {len(tokeep)}")

## Branch the genealogies
new_ped = gen.branching(ped, pro=tokeep)
new_asc = gen.genout(new_ped)

## Save the branch
out_file = sys.argv[2]
new_asc.to_csv(
        out_file,
        sep='\t',
        index=False
)
