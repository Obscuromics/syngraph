suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
# require('ggtree')

parser <- ArgumentParser()
parser$add_argument("-i", "--input",  
    dest="i", help="Input tree to plot, use syngraph output so the node labels correspond to the inference (with .newick.txt suffix).")
parser$add_argument("-o", "--output",  
    dest="o", help="Name of the output pdf tree (including .pdf suffix)")
parser$add_argument("-c", "--clusters",  
    dest="c", help="Name of (*.clusters.tsv)", default = "")

args <- parser$parse_args()

# clusers_file <- 'data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.clusters.tsv'
clusers_file <- args$c
# tree_file <- "data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.newick.txt"
tree_file <- args$i

tree <- read.tree(tree_file)

if (clusers_file != ""){
    clusters <- read.table(clusers_file, col.names = c('clust', 'node', 'LGs'))
    only_those_in_the_tree <- clusters[clusters[, 'node'] %in% tree$node.label, ]
    row.names(only_those_in_the_tree) <- only_those_in_the_tree[, 'node']
    tree$node.label <- paste(tree$node.label, ":", only_those_in_the_tree[tree$node.label, 'LGs'])
    tree$node.label[1] <- ""
}

pdf(args$o,height = 60, width = 20)
    plot(tree,  show.node.label = TRUE)
dev.off()
