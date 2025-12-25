from bx.align import axt
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import gzip
import numpy as np

# edit this according to where you keep the ucsc alignment
AXT_CHR_PATH = lambda chr_id : f"../ucsc-rat-mus/axtNet/chr{chr_id}.mm39.rn7.net.axt.gz"

def get_rat_seq(coords, chr_id):
    """
    Get the sequence of rat characters corresponding
    to mouse in the whole genome alignment and the mask
    marking which characters were not aligned.
    """
    handle = gzip.open(AXT_CHR_PATH(chr_id), 'rt')
    axt_reader = axt.Reader(handle)

    aligned_mask = np.zeros(coords.size, dtype=bool)
    # mus_seq = np.full(len(coords), '?')
    rat_seq = np.full(len(coords), '?')

    curr_aln = next(axt_reader)
    mus, rat = curr_aln.components
    st, ed = mus.start, mus.end
    nonindel = np.array([x for x in mus.text]) != '!'
    seq_indexer = np.arange(len(mus.text))

    N = len(coords)
    ix = 0
    while ix < N:
        coord = coords[ix]
        if coord < st: # coordinate appears before start, move to next coordinate
            ix += 1
        else: # coordinate appears after start
            if coord < ed: # coordinate appears within range
                relative_pos = coord - mus.start
                true_rel_pos = seq_indexer[nonindel][relative_pos]
                # mus_seq[ix] = mus.text[true_rel_pos]
                rat_seq[ix] = rat.text[true_rel_pos]
                aligned_mask[ix] = True
                ix += 1
            else: # coordinate appears after end, move to next alignment
                curr_aln = next(axt_reader)
                if not curr_aln:
                    break
                mus, rat = curr_aln.components
                st, ed = mus.start, mus.end
                nonindel = np.array([x for x in mus.text]) != '-'
                seq_indexer = np.arange(len(mus.text))

    # Sanity check with mouse sequence:
    # return ''.join(mus_seq), aligned_mask
    
    return rat_seq, aligned_mask

def add_rat_seq(msa, coords, chr_id, filter_gaps=False):
    """
    Return an MSA object with rat's sequence added as an outgroup.
    Unaligned sites (marked by ?'s) are omitted in the combined alignment.
    """
    rat_seq, aligned_mask = get_rat_seq(coords, chr_id)
    if filter_gaps:
        aligned_mask = (rat_seq != '-') & aligned_mask
    rat_seq = Seq(''.join(rat_seq[aligned_mask]).upper())
    rat_record = SeqRecord(rat_seq, id="rat", description="")
    result = MultipleSeqAlignment([rat_record])
    
    for seq_record in msa:
        masked_seq = Seq("".join(np.array(seq_record.seq)[aligned_mask]))
        masked_record = SeqRecord(masked_seq, id=seq_record.id, description="")
        result.append(masked_record)

    return result, coords[aligned_mask]


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("coordinates")
    parser.add_argument("chr_id")
    parser.add_argument("msa_path")
    
    args = parser.parse_args()

    # AXT coordinates are 0-indexed
    # VCF coordinates are 1-indexed
    coords = np.load(args.coordinates) - 1 # If the coordinates are 1-indexed, subtract by 1
    chr_id = args.chr_id
    msa = AlignIO.read(args.msa_path, "fasta")
    
    msa_with_rat, masked_coords = add_rat_seq(msa, coords, chr_id, filter_gaps=True)
    AlignIO.write(msa_with_rat, f"{args.msa_path.split('.fasta')[0]}_rat.fasta", "fasta")

    # save the filtered coordinates in the original base index
    np.save(f"{args.coordinates.split('.npy')[0]}_rat.npy", masked_coords+1)

if __name__ == "__main__":
    main()

