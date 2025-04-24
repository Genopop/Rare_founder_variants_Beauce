## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import sys
from sys import argv

print("Describe the final regional groups")
## Get the pedigree informations
ped_file = sys.argv[1]
ped = gen.genealogy(ped_file)
## Number of probands
pro = gen.pro(ped)
print(f"Number of probands: {len(pro)}")
## Number of ancestors
ancestors = gen.ancestor(ped, pro)
print(f"Number of ancestors: {len(ancestors)}")
## Number of founders
founders = gen.founder(ped)
print(f"Number of founder: {len(founders)}")
## Maximum depth
depth = gen.depth(ped)
print(f"Maximum depth: {depth}")
