#!/usr/bin/env python3

import sys
from os import listdir
from os.path import isfile, join
from Bio import Phylo

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            next
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

treefile = './data/results/infer_alg_diptera/iqtree/diptera.supermatrix.phy.treefile'
exclusion_dir = './data/results/infer_alg_diptera/excluded_syngraph_busco_tables/'
outgroup_f = None # './data/results/infer_alg_diptera/outgroup_clade.txt'

exclusion_files = [f for f in listdir(exclusion_dir) if isfile(join(exclusion_dir, f))]
whole_exclusions = [f.split('.')[0] for f in exclusion_files if 'syngraph.buscos.tsv' in f]

tree = Phylo.read(treefile, "newick")
tree_names = tabulate_names(tree)
leaves = tree.get_terminals()
excluded_leaves = [leaf for leaf in leaves if leaf.name in whole_exclusions]

for leaf in excluded_leaves:
	tree.prune(leaf)

if outgroup_f:
	with open(outgroup_f, 'r') as f:
		outgroup_taxa = f.readlines()

	outgroup = [{"name": taxon} for taxon in outgroup_taxa]

	if tree.is_monophyletic(outgroup):
		tree.root_with_outgroup(*outgroup)
		Phylo.write(tree, './data/results/infer_alg_diptera/diptera.pruned.rerooted.newick','newick')
	else:
    		raise ValueError("outgroup is paraphyletic")

else:
	Phylo.write(tree, './data/results/infer_alg_diptera/diptera.pruned.newick','newick')
