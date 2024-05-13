#!/usr/bin/env python3
import argparse
import numpy as np
import ete3

def scale_tree(tree, x, inplace=False):
    if inplace:
        result_tree = tree
    else:
        result_tree = tree.copy()
    for e in result_tree.traverse():
        e.dist *= x
    return result_tree

def nonultrametric_rescaling(tree, c, inplace=False):
    if inplace:
        result_tree = tree
    else:
        result_tree = tree.copy()
    logc = np.log(c)
    for e,x in zip(result_tree.traverse(), np.exp(np.random.uniform(-logc, logc, 2*len(result_tree)-1))):
        e.dist *= x
    return result_tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scale branches")

    parser.add_argument(
        "trees", type=ete3.Tree, nargs="*",
        help="Tree in the newick format"
    )
    parser.add_argument(
        "-c", "--deviation_factor", type=float,
        help="Deviation factor for introducing nonultrametricity."
    )
    parser.add_argument(
        "-s", "--seed", type=int,
        help="Random number generator seed"
    )
    parser.add_argument(
        "-H", "--height", type=float,
        help="Scale tree height to specified value"
    )
    parser.add_argument(
        "-m", "--multiply", type=float,
        help="Multiply branch lengths by specified value"
    )

    args = parser.parse_args()

    c = args.deviation_factor
    if args.seed:
        np.random.seed(args.seed)

    for tree in args.trees:
        if args.height:
            if args.multiply:
                print("Cannot use both -h/--height and -m/--multiplty")
                exit(1)
            multiplier = args.height / tree.get_farthest_leaf()[1]
            scale_tree(tree, multiplier, inplace=True)
        elif args.multiply:
            scale_tree(tree, args.multiply, inplace=True)
        if args.deviation_factor:
            nonultrametric_rescaling(tree, args.deviation_factor, inplace=True)
        print(tree.write(format=1))
