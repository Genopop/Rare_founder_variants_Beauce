## Import packages
import geneo as gen
import pandas as pd
import numpy as np
import pickle
import sys
from sys import argv

pickle_file = sys.argv[1]
with open(pickle_file, "rb") as f:
    phi = pickle.load(f)

ci = gen.phiCI(phi)
CI_file = sys.argv[2]
ci.to_csv(CI_file)
