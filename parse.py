from tqdm import tqdm
import argparse
import re

# parser = argparse.ArgumentParser()
# parser.add_argument('xmfa', type=str)
# parser.add_argument('fna', type=str)
# args = parser.parse_args()

xmfa_path = "./parsnp.xmfa"
ref_path = "./"
fna_path = "./test/"

with open(xmfa_path) as xmfa:
    line = xmfa.readline()
    line = xmfa.readline()
    seqNum = int(line.split()[1])
    seqs = {}
    for i in tqdm(range(seqNum)):
        line = xmfa.readline()
        line = xmfa.readline()
        seqName = line.split()[1]
        if seqName[-4:] == ".ref":
            seqName = seqName[:-4]
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
            while line and (initial:=line[0]) != '>':
                if initial == '=':
                    pbar.update(1)
                line = xmfa.readline()

# parsing xmfa done.

for seq, coords in tqdm(seqVerify.items()):
    if seq == 1:
        path = ref_path + seqs[seq]
    else:
        path = fna_path + seqs[seq]
    coords.sort(key=lambda x:x[0])
    with open(path) as fna:
        it_coords = iter(coords)
        counter = 0
        target, correct = next(it_coords, (None, None))
        line = fna.readline()
        while line and target:
            if counter >= target:
                # if we hit over the target, get the relative position and start 
                # comparing fna (`compare`) to xmfa (`correct`)
                pos = target - (counter - (len(line) - 1))
                compare = line[pos:].strip()
                if (length:=len(compare)) > 5:
                    cutoff = min(20, length)
                    if correct[:cutoff].lower() != compare[:cutoff].lower():
                        print(counter)
                        # print(length)
                        print(compare[:cutoff].lower())
                        print(correct[:cutoff].lower())
                        print(target)
                        print("----")
                        break
            
                target, correct = next(it_coords, (None, None))
                # updating the target
            else:
                line = fna.readline()
                if line[0] != '>':
                    counter += len(line.strip())
                    # updating the counter for every non-header line