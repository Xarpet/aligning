from tqdm import tqdm
import argparse
import re

# parser = argparse.ArgumentParser()
# parser.add_argument('xmfa', type=str)
# parser.add_argument('fna', type=str)
# args = parser.parse_args()

xmfa_path = "parsnp.xmfa"
fna_path = "./"

with open(xmfa_path) as xmfa:
    line = xmfa.readline()
    line = xmfa.readline()
    seqNum = int(line.split()[1])
    seqs = {}
    for i in tqdm(range(seqNum)):
        line = xmfa.readline()
        line = xmfa.readline()
        seqName = line.split()[1]
        seqs[i+1] = seqName
        line = xmfa.readline()
        line = xmfa.readline()

    print(seqs)

    line = xmfa.readline()
    intervalCount = int(line.split()[1])
    print(intervalCount)

    # Header parsing over

    seqVerify = {}
    for seq in range(seqNum):
        seqVerify[seq+1] = []

    line = xmfa.readline()

    with tqdm(total=intervalCount) as pbar:
        while line:
            alignment = re.split("-|:p| cluster| s|:|\s", line[1:])
            # Here the alignments are in order: 
            # [seqeunce number, starting coord, end coord, + or -, cluster number, contig number, coord in contig]
            line = xmfa.readline()
            if alignment[3] == "+":
                seqVerify[int(alignment[0])].append((int(alignment[1]),line[:20]))
            # notice that here only the first 20 are taken
            # get to next alignment header
            while(line and (initial:=line[0]) != '>'):
                if initial == '=':
                    pbar.update(1)
                line = xmfa.readline()

# parsing xmfa done.

for seq in seqVerify:
    path = seqs[seq]
    with open(fna_path+path) as fna:
        print("1")
