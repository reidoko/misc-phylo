#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from pathlib import Path
import numpy as np
from collections import Counter
from itertools import combinations

def char_freqs(aln):
    result = ""
    counts = Counter(aln.flat)
    if '-' in counts:
        counts.pop('-')
    for ch in sorted(counts.keys()):
        result = result + f"\n  {ch.upper()}: {(100*counts[ch]/counts.total()):.1f} "
    return result[:-1]

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
    description="Calculate basic summary statistics about an alignment (length, gappiness, ANHD)"
)
parser.add_argument(
    "alignment", type=read_alignment, help="Path to alignment"
)

args = parser.parse_args()

aln = np.array(args.alignment)
nseqs, aln_length = aln.shape
anhds = []
for x,y in combinations(aln, 2):
    anhds.append(((x!=y)&((x!='-')|(y!='-'))).sum() / aln_length)
gappis = []
for x in aln:
    gappis.append((x=='-').sum() / aln_length)

anhds = np.array(anhds)
gappis = np.array(gappis)

print(f"""Length: {aln_length}
ANHD: {anhds.mean():.3f}
Gappiness: {gappis.mean():.3f}
Character frequencies: {char_freqs(aln)}""")
