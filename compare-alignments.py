#!/usr/bin/env python3
import argparse
from Bio import AlignIO
from pathlib import Path
import numpy as np
from itertools import combinations
from collections import Counter

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


def SP_FN_FP_jacc(true_alignment, inferred_alignment):
    # Assumes sequences appear in the same order, i.e. taxon1 is the first taxon for both true and inferred
    true_al = true_alignment != '-'
    infd_al = inferred_alignment != '-'
    n = true_al.shape[0]
    taxa = list(range(n))
    
    total_true_homologies = 0
    total_est_homologies = 0
    total_false_positive = 0
    total_false_negative = 0

    jaccard_numer = 0
    jaccard_denom = 0
    
    for ix,jx in combinations(taxa, 2):
        taxa_mask = np.isin(taxa, (ix,jx))
        t_homologies = (true_al[taxa_mask].cumsum(axis=1)-1).T[true_al[ix] & true_al[jx]]
        total_true_homologies += t_homologies.shape[0]
        true_homologies = frozenset(map(lambda x: (x[0], x[1]), t_homologies))
         
        i_homologies = (infd_al[taxa_mask].cumsum(axis=1)-1).T[infd_al[ix] & infd_al[jx]]
        total_est_homologies += i_homologies.shape[0]
        inferred_homologies = frozenset(map(lambda x: (x[0], x[1]), i_homologies))
        
        total_false_positive += len(inferred_homologies - true_homologies)
        total_false_negative += len(true_homologies - inferred_homologies)
        jaccard_numer += len(inferred_homologies.symmetric_difference(true_homologies))
        jaccard_denom += len(inferred_homologies.union(true_homologies))
    return total_false_negative / total_true_homologies, total_false_positive / total_est_homologies, jaccard_numer / jaccard_denom

def bigfoot_encode(seq):
    # Based on encoding in BigFoot supplementary materials
    # https://doi.org/10.1186/1471-2148-9-217
    not_gaps = (np.array(seq) != '-')
    sumseq = not_gaps.cumsum()*2-1
    sumseq[~not_gaps] += 1
    return sumseq

def bigfoot_cols(A):
    return np.array([bigfoot_encode(seq.seq) for seq in A])

def TC_score(al1, al2):
    # Following the definition in https://doi.org/10.2478/v10006-009-0054-y
    # 
    # 
    # gaps all use the same symbol (0),
    # otherwise, the second and third columns are not counted.
    # AAAA = 1357 -> 1357
    # A--- = 1222 -> 1000
    # and
    # AAAA = 1357 -> 1357
    # ---A = 0001 -> 0001
    bf1 = bigfoot_cols(al1)
    bf2 = bigfoot_cols(al2)
    bf1[(bf1+1) % 2] = 0
    bf2[(bf2+1) % 2] = 0
    cols1 = Counter(tuple(bf1[:,col]) for col in range(bf1.shape[1]))
    cols2 = Counter(tuple(bf2[:,col]) for col in range(bf2.shape[1]))
    
    return (cols1 & cols2).total()

def summary_stats(aln):
    nseqs, aln_length = aln.shape
    anhds = []
    for x,y in combinations(aln, 2):
        anhds.append(((x!=y)&((x!='-')|(y!='-'))).sum() / aln_length)
    gappis = []
    for x in aln:
        gappis.append((x=='-').sum() / aln_length)

    anhds = np.array(anhds)
    gappis = np.array(gappis)

    print(f"""    Length: {aln_length}
    ANHD: {anhds.mean():.3f}
    Gappiness: {gappis.mean():.3f}""")


parser = argparse.ArgumentParser(
    description="Get SP-FN, SP-FP, TC, and jaccard of two alignments."
)
parser.add_argument(
    "alignment1", type=read_alignment, help="path to reference alignment")
parser.add_argument(
    "alignment2", type=read_alignment, help="path to other alignment")


args = parser.parse_args()

aln1 = np.array(args.alignment1)
aln2 = np.array(args.alignment2)


print("Reference alignment:")
summary_stats(aln1)
print("Estimated alignment:")
summary_stats(aln2)

spfn, spfp, jacc = SP_FN_FP_jacc(aln1, aln2)
tc_score = TC_score(args.alignment1, args.alignment2)
print(f"""
SP-FP: {spfp:.3f}
SP-FN: {spfn:.3f}
TC score: {tc_score}
Jaccard distance: {jacc:.3f}""")

