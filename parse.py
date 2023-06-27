from functools import reduce
import operator
import os
import re
from tqdm import tqdm
import argparse
from datetime import datetime
from Bio import SeqIO
from Bio import AlignIO

# python parse.py parsnp.xmfa ./ ./test/ -m

# I made them positional arguments so that there's no need of dashes
# also I feel like we could just use the file names in xmfa? that way we
# only need to provide path to all the fna files; not sure
parser = argparse.ArgumentParser()
parser.add_argument('xmfa', type=str, help="the path to the xmfa file you are trying to verify")
parser.add_argument('ref', type=str, help="the path to the referrence fna file")
parser.add_argument('fna', type=str, help="the path to all the fna files")
parser.add_argument('-m', '--maf', action='store_true', help="exports a .maf file translated from the xmfa file")
args = parser.parse_args()

xmfa_path = args.xmfa
ref_path = args.ref
fna_path = args.fna
maf_flag = args.maf

def compare_with_dashes(str1, str2):
    # ignores the dashes when comparing
    if str1 == str2:
        return True
    if len(str1) != len(str2):
        return False
    else:
        return all(c1 == c2 for (c1, c2) in filter(lambda pair: '-' not in pair, zip(str1, str2)))

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
                # here only forward alignments are used
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
        dna = str(reduce(operator.add, [record.seq for record in SeqIO.parse(path, "fasta")]))
        # flatmapping all sequence together
        for target, correct in coords:
            length = len(correct)
            if not compare_with_dashes(compare:=dna[target:target+length].lower(), correct.lower()):
                f.write("sequence: " + str(seq) + "\n")
                f.write("file name: " + seqs[seq] + "\n")
                f.write("position: " + str(target) + "\n")
                if (actual_pos:=dna.lower().find(correct)) != (-1):
                    f.write("actual position: " + str(actual_pos) + '\n')
                f.write("fna: " + compare + "\n")
                f.write("xmfa: " + correct.lower() + "\n")
                f.write("----" + "\n")


if maf_flag:
    data = []
    with open(xmfa_path, "r+") as xmfa:
        for line in xmfa:
            if line[0] == '>' and line[1] != " ":
                data.append("> " + line[1:])
            else:
                data.append(line)

    with open("temp_xmfa", 'x') as tmp:
        for d in data:
            tmp.write(d)

    alignments = AlignIO.parse("temp_xmfa", 'mauve')

    with open(current_time+".maf", "x") as file:
        maf = AlignIO.MafIO.MafWriter(file)
        for align in alignments:
            maf.write_alignment(align)

    os.remove("temp_xmfa")


# legacy code
"""
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
"""
