#!/usr/bin/env python3
import argparse
from ete3 import Tree
from itertools import combinations
from collections import Counter, defaultdict
from pathlib import Path

def get_leafset(node):
    if node.is_leaf():
        return frozenset({node.name})
    else:
        return frozenset({x.name for x in node.get_leaves()})
    
def get_bipartition(node, full_leafset):
    left = get_leafset(node)
    right = full_leafset - left
    return frozenset({left, right})

def is_trivial(bip):
    left, right = bip
    return len(left) <= 1 or len(right) <= 1

def get_nontrivial_bipartitions(tree):
    biparts = set()
    all_leaves = frozenset([x.name for x in tree.get_leaves()])
    for node in tree.traverse():
        if not node.is_leaf():
            x = get_bipartition(node, all_leaves)
            if not is_trivial(x):
                biparts.add(x)

    return frozenset(biparts)

def mean_discordance(trees):
    # there's probably a faster way that just looks at the counts of the bipartitions
    tree_counts = Counter(map(get_nontrivial_bipartitions, trees))
    n_trees = len(trees)
    discordance = 0
    for (b1,c1),(b2,c2) in combinations(tree_counts.items(), 2):
        discordance += c1*c2*(len(b1.symmetric_difference(b2))/(len(b1) + len(b2))) 
    return discordance / (n_trees * (n_trees-1)//2)

def newick_or_file(inp):
    input_path = Path(inp)
    if input_path.exists():
        return [Tree(x) for x in open(input_path) if len(x) > 0]
    else:
        return [Tree(newick=inp)]
    
def flatten(list):
    return [x if not isinstance(sublist, Tree) else sublist for sublist in list for x in sublist ]
parser = argparse.ArgumentParser(description="Calculate NRF between all combinations of trees in a file")
parser.add_argument("trees", nargs="*", type=newick_or_file, help="Trees to use")

args = parser.parse_args()
trees = flatten(args.trees)
print(mean_discordance(trees))
