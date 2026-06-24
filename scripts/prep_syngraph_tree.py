#!/usr/bin/env python3

import sys
from os import listdir
from os.path import isfile, join
from Bio import Phylo
import argparse

parser = argparse.ArgumentParser(description="Prepare a syngraph tree by pruning the tree to keep only the taxa present in the syngraph table, and optionally rooting the tree with an outgroup")
parser.add_argument("-intree", required=True, help="Path to the input tree file (Newick format)")
parser.add_argument("-incl_dir", required=True, help="Path to the directory with per species busco tables formated for syngraph inference.")
parser.add_argument("-outtree", required=True, help="Path to the output tree file")
parser.add_argument("-root", help="Tip name of the root (default: don't root)", default='')

args = parser.parse_args()

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            next
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

treefile = args.intree # './data/results/infer_alg_diptera/iqtree/diptera.supermatrix.phy.treefile'
outtree = args.outtree
# exclusion_dir = sys.argv[2] # './data/results/infer_alg_diptera/excluded_syngraph_busco_tables/'
inclusion_dir = args.incl_dir # './data/results/infer_alg_diptera/excluded_syngraph_busco_tables/'

outgroup_f = args.root # './data/results/infer_alg_diptera/outgroup_clade.txt'

inclusion_files = [f for f in listdir(inclusion_dir) if isfile(join(inclusion_dir, f))]
whole_inclusions = [f.split('.')[0] for f in inclusion_files if 'syngraph.buscos.tsv' in f]

tree = Phylo.read(treefile, "newick")
tree_names = tabulate_names(tree)
leaves = tree.get_terminals()
excluded_leaves = [leaf for leaf in leaves if leaf.name not in whole_inclusions]

for leaf in excluded_leaves:
	tree.prune(leaf)

# print(outgroup_f)
# if outgroup_f:
# 	# if isfile(outgroup_f):
# 	# 	with open(outgroup_f, 'r') as f:
# 	# 		outgroup_taxa = f.readlines()
# 	# else:
# 	# 	outgroup_taxa = [outgroup_f]
# 	# outgroup = [{"name": taxon.strip()} for taxon in outgroup_taxa]
# 	# if tree.is_monophyletic(outgroup):
# 	tree.root_with_outgroup(outgroup_f)

Phylo.write(tree, outtree,'newick')
