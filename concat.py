#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from functools import reduce
from pathlib import Path
from tqdm import tqdm

def read_alignment(a_path, form=None):
    a_path = Path(a_path)
    suffix = a_path.suffix.lower()
    valid_alignment_formats = {
        "clustal",
        "emboss",
        "fasta",
        "fasta-m10",
        "ig",
        "maf",
        "mauve",
        "msf",
        "nexus",
        "phylip",
        "phylip-sequential",
        "phylip-relaxed",
        "stockholm"
    }
    
    if form:
        if form not in valid_alignment_formats:
            raise ValueError(f"{form} not a valid alignment format specification, must be one of the following: [ {','.join(valid_argument_formats)} ]")
    else:
        if not suffix:
            raise NotImplementedError("Not able to guess alignment type without a suffix, provide a format.")
        elif suffix in {".ph", ".phy", ".phylip"}:
            form="phylip-relaxed"
        elif suffix in {".fa", ".fna", ".ffn", ".faa", ".frn", ".fasta"}:
            form="fasta"
        elif suffix in {".nx", ".nex", ".nxs"}:
            form="nexus"
        else:
            raise NotImplementedError(f"Not able to resolve alignment suffix {suffix} yet")
    
    # if a_path.exists():
    al = AlignIO.read(a_path, format=form)
    al.sort()
    return al

parser = argparse.ArgumentParser(
    description="Concatenate multiple sequence alignments with identical taxon sets"
)
parser.add_argument(
    "-a", "--alignments", type=Path, nargs="*",
    help="Paths to alignments to concatenate"
)
parser.add_argument(
    "-o", "--output", type=Path, required=True,
    help="Output path"
)
parser.add_argument(
    "-f", "--format", type=str,
    help="Output alignment format", default="fasta"
)
parser.add_argument(
    "-p", "--partitions", type=str, default=None,
    help="Partition file path for RAxML-NG or NetRax"
)
parser.add_argument(
    "-m", "--model", type=str, default="GTR+G",
    help="Partition file model type"
)

args = parser.parse_args()
alignments = map(read_alignment, args.alignments)
# concatenated = reduce(lambda x,y: x+y, alignments)
concatenated = next(alignments)
ix = 0
start=1
partitions = [f"{args.model}, {args.alignments[ix].stem} = {start}-{concatenated.get_alignment_length()}"]
# print(args.alignments[ix].stem)
for msa in tqdm(alignments, total=len(alignments)):
    ix += 1
    start = concatenated.get_alignment_length()
    concatenated += msa
    partitions += [f"{args.model}, {args.alignments[ix].stem} = {start+1}-{concatenated.get_alignment_length()}"]
AlignIO.write(concatenated, args.output, args.format)
if args.partitions:
    with open(args.partitions, 'w') as handle:
        handle.write('\n'.join(partitions)+'\n')
