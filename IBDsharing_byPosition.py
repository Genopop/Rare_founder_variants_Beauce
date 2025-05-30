import numpy as np
import sys
from sys import argv

# Arguments: path and file names
# datafile: .merged.ibd files (Carriers only)
# mychr : chromosome 
# bim_file : PLINK bim file
datafilename = sys.argv[1]
mychr = int(sys.argv[2])
outputfilename = datafilename + '.sharing.by.pos'
bim_file = sys.argv[3]

# Load IBD sharing file
def loadfile(datafile):
    SNP1 = np.genfromtxt(fname=datafile, dtype=int, usecols=5)
    SNP2 = np.genfromtxt(fname=datafile, dtype=int, usecols=6)
    return [SNP1, SNP2]

# Load BIM file for map
def loadmapfile(mapfile):
    positions = []
    with open(mapfile, 'r') as f:
        for line in f:
            ichr = int(line.rstrip('\n\r').split()[0])
            if ichr == mychr:
                pos = int(line.rstrip('\n\r').split()[3])
                positions.append(pos)
    return positions

# Write file
def writefile(outputfilename, dict, pos):
    with open(outputfilename, "w") as outputfile:
        for i in range(len(pos)):
            outputfile.write(f"{mychr}\t{pos[i]}\t{dict[pos[i]]}\n")

# Main
print('loading data file...')
data = loadfile(datafilename)

positions = loadmapfile(bim_file)
print(len(positions))
print(positions[1:10])
print("Positions done...")

pos_dict = dict(zip(positions, [0] * len(positions)))
print("Dict done...looping through IBD segments...")

for i in range(len(data[0])):
    try:
        start_index = positions.index(data[0][i])
        end_index = positions.index(data[1][i]) + 1
        print(f"Processing segment: {data[0][i]} to {data[1][i]} (indexes {start_index} to {end_index})")
        for j in range(start_index, end_index):
            pos_dict[positions[j]] += 1
            print(f"Incremented position {positions[j]}: count now {pos_dict[positions[j]]}")
    except ValueError as e:
        print(f"Error: {e}. Skipping segment [{data[0][i]}, {data[1][i]}]")

writefile(outputfilename, pos_dict, positions)
