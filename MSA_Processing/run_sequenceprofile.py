import Bio
import Bio.PDB

import argparse

parser = argparse.ArgumentParser(
    description="Transform the MSA information to the sequence profile distribution"
)
parser.add_argument("-i", "--input", type=str, required=True, help="Input file (aln format)")
parser.add_argument(
    "-s",
    "--scale",
    type=float,
    default=1,
    help="Scaling factor. Utilizing smaller values leads to a sharper distribution",
)
parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
args = parser.parse_args()

AA = {
    "A": 0,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 5,
    "H": 6,
    "I": 7,
    "K": 8,
    "L": 9,
    "M": 10,
    "N": 11,
    "P": 12,
    "Q": 13,
    "R": 14,
    "S": 15,
    "T": 16,
    "V": 17,
    "W": 18,
    "Y": 19,
    "-": 20,
    "U": 20,
    "X": 20,
    "B": 20,
    "Z": 20,
}

f = open(args.input, "r")
rows = f.readlines()
f.close()

sequence_profile = []
for row in rows:
    seqs = []
    for i in row.strip():
        seqs.append(AA[i])
    sequence_profile.append(seqs)

import numpy as np
import torch
import torch.nn.functional as F

data = torch.tensor(sequence_profile).transpose(0, 1)
distri = []
for i in range(data.size(0)):
    hist = torch.histc(data[i].float(), bins=21, min=0, max=20)
    prob_dist = hist / hist.sum()
    distri.append(prob_dist)

data = torch.stack(distri)
print("Sequence profile:", data.size())
if args.scale == 1:
    torch.save(data, args.output)
else:
    scaled_data = data / args.scale
    torch.save(F.softmax(scaled_data, dim=1), args.output)
