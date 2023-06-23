from tqdm import tqdm
import argparse
import re
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument('xmfa', type=str, help="the path to the xmfa file you are trying to verify")
parser.add_argument('ref', type=str, help="the path to the referrence fna file")
parser.add_argument('fna', type=str, help="the path to all the fna files")
args = parser.parse_args()

xmfa_path = args.xmfa
ref_path = args.ref
fna_path = args.fna

def compare_with_dashes(str1, str2):
    # ignores the dashes when comparing
    if str1 == str2:
        return True
    if len(str1) != len(str2):
        return False
    for index in range(len(str1)):
        if (char1:=str1[index]) != (char2:=str2[index]):
            if char1 != '-' and char2 != '-':
                return False
    return True

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

    line = xmfa.readline()
    intervalCount = int(line.split()[1])

    # Header parsing over

    seqVerify = {}
    for seq in range(seqNum):
        seqVerify[seq+1] = []

    line = xmfa.readline()

    with tqdm(total=intervalCount) as pbar:
        while line:
            alignment = re.split("-|:p| cluster| s|:|\s", line[1:])
            # Here the alignments are in order:
            # [seqeunce number, starting coord, end coord, ...
            # + or -, cluster number, contig number, coord in contig]
            line = xmfa.readline()
            if alignment[3] == "+":
                seqVerify[int(alignment[0])].append(
                    (int(alignment[1]), line[:20]))
            # notice that here only the first 20 are taken
            # get to next alignment header
            while line and (initial := line[0]) != '>':
                if initial == '=':
                    pbar.update(1)
                line = xmfa.readline()

# parsing xmfa done.

now = datetime.now()

current_time = now.strftime("%Y-%m-%d-%H%M%S")

with open(current_time+".txt", "x") as f:
    for seq, coords in tqdm(seqVerify.items()):
        if seq == 1:
            path = ref_path + seqs[seq]
        else:
            path = fna_path + seqs[seq]
        coords.sort(key=lambda x: x[0])
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
                    if (length := len(compare)) > 5:
                        cutoff = min(20, length)
                        if not compare_with_dashes(correct[:cutoff].lower(), compare[:cutoff].lower()):
                            # print(counter)
                            # print(length)
                            f.write("sequence: " + str(seq) + "\n")
                            f.write("file name: " + seqs[seq] + "\n")
                            f.write("position: " + str(target) + "\n")
                            f.write("fna: " + compare[:cutoff].lower() + "\n")
                            f.write("xmfa: " + correct[:cutoff].lower() + "\n")
                            f.write("----" + "\n")

                    target, correct = next(it_coords, (None, None))
                    # updating the target
                else:
                    line = fna.readline()
                    if line[0] != '>':
                        counter += len(line.strip())
                        # updating the counter for every non-header line
