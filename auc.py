#!/usr/bin/env python3

import ete3
import argparse

from pathlib import Path
from functools import reduce
from collections import Counter

from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay

import matplotlib.pyplot as plt

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

    return biparts

def get_raxml_nontrivial_bips_and_support(tree):
    biparts = set()
    bs = []
    all_leaves = frozenset([x.name for x in tree.get_leaves()])
    for node in tree.traverse():
        if not node.is_leaf():
            x = get_bipartition(node, all_leaves)
            trivial = is_trivial(x)
            if x not in biparts and not trivial:
                bs.append((x, node.support))
            biparts.add(x)

    return bs

def branch_support(tree, reests):
    """
    Args:
        tree: ete3 Tree
        reests: list of ete3 Trees
    
    Returns:
        List of pairs with bipartition and the corresponding support value (i.e # bipartitions present in re-estimated trees / # of re-estimated trees)
    """
    ref_bips = get_nontrivial_bipartitions(tree)
    supp_bips = reduce(lambda x,y : list(x)+list(y), map(get_nontrivial_bipartitions, reests))
    bip_counts = Counter(supp_bips)
    n = len(reests)
    
    return [(bip, bip_counts[bip]/n) for bip in ref_bips]

def plot_curves(truth, pred):
    
    precision, recall, thresholds = precision_recall_curve(truth, pred)
    fpr, tpr, thresholds = roc_curve(truth, pred)
    roc_auc = roc_auc_score(truth, pred)
    
    fig, axs = plt.subplots(1,2)

    axs[0].set_title("Precision-Recall Curve")
    pr_display = PrecisionRecallDisplay(precision=precision, recall=recall)
    pr_display.plot(ax=axs[0])
    
    roc_display = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
    roc_display.plot(ax=axs[1])
    axs[1].set_title("ROC Curve")

    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', type=ete3.Tree, nargs='+', help="Model tree in newick format")
    parser.add_argument('-e', '--estimated', type=ete3.Tree, nargs='+', help="Estimated tree in newick format")
    parser.add_argument('-r', '--reestimate', type=Path, nargs='+', help="Path to file containing re-estimated trees")
    parser.add_argument('-p', '--plot', action='store_true', help='Set to plot PR and ROC curves')

    args = parser.parse_args()

    if not (args.model and args.estimated and args.reestimate):
        parser.print_help()
        exit(1)

    reests = []
    for reest_path in args.reestimate:
        with open(reest_path, 'r') as handle:
            reest_trees = list(map(ete3.Tree, filter(lambda x: len(x.strip())>0, handle.readlines())))
        reests.append(reest_trees)
    if len(args.model) != len(args.estimated) or len(args.model) != len(reests):
        print(f"Error: Not the same number of model trees, estimated trees, and re-estimated tree files! Only found {len(args.model)} model trees, {len(args.estimated)} estimated trees, and {len(reests)} re-estimated tree files.")
        exit(1)
    
    truth = []
    pred = []

    for mtree, etree, reest in zip(args.model, args.estimated, reests):
        supports = branch_support(etree,reest)
        true_bips = get_nontrivial_bipartitions(mtree)
        for bip, sup in supports:
            truth.append(bip in true_bips)
            pred.append(sup)

    # I don't remember why I didn't use this -- it should be the same?
    # fpr, tpr, thresholds = roc_curve(truth, pred)
    # roc_auc = auc(fpr, tpr)
    precision, recall, thresholds = precision_recall_curve(truth, pred)
    pr_auc = auc(recall, precision)
    roc_auc = roc_auc_score(truth, pred)

    if args.plot:
        plot_curves(truth,pred)

    print(f"PR-AUC:  {pr_auc}\nROC-AUC: {roc_auc}")

if __name__ == "__main__":
    main()
