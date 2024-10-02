# Miscellaneous scripts, notes, etc for phylogenetics

## Scripts
General requirements to run any of the python scripts are in requirements.txt

**auc.py**  
Calculate AUC for precision-recall and receiver operating characteristic curves
for tree topology support estimates.

**concat.py**  
Concatenate multiple sequence alignments for use in other analyses (e.g. RAxML concatenated).

**pairwise-nrf.py**  
Calculate the mean robinson-foulds distance between each pair of trees given as input.

**scaletree.py**  
Scale branch-lengths either by a constant factor or to achieve a desired tree height 
(maximum root-to-tip distance). Can also introduce ultrametricity by the procedure 
described in [Moret et al., 2002](https://doi.org/10.1007/3-540-45784-4_26)
