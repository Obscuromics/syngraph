suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
# suppressPackageStartupMessages(library(ggplot2))

# parser <- ArgumentParser()

# parser$add_argument("-a", "--alg", 
#     help="The name of the ALG definitions to be used.")
# parser$add_argument("-o", "--output",  
#     dest="o", help="Base of the output name (.png will be attached)")

# args <- parser$parse_args()
source('scripts/20250620_colour_pal.R')

# these need to be made into arguments in the end
treefile <- 'data/syngraph/diptera.no_plecia.mindist.m165.newick.txt'
rearrangement_file <- 'data/syngraph/diptera.no_plecia.mindist.m165.rearrangements.tsv'
busco_asn_file <- 'tables/ALGs_syngraph_diptera.tsv' # args$a

output <- 'figures/syngraph_tree_of_changes.pdf'
# simplified; if 0 means not simplified, if a number it is the number of pieces it gets simplfied to
simplified <- 0
just_internal <- TRUE

directory_with_nodes <- paste0("data/syngraph/node_assignments/")
all_node_files <- dir(directory_with_nodes)

ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'chrQ'))

########## Loading states at all the nodes
# testing: input_file <- paste0("data/syngraph/node_assignments/", all_node_files[1])
load_node <- function(node_file){
    node_lgs <- read.table(node_file, header = T, col.names = c('busco', 'chrNode'))
    # needs ALG_data
    node_lgs_by_algs <- merge(ALG_data, node_lgs)
    dat <- data.frame(table(node_lgs_by_algs$chrQ,node_lgs_by_algs$chrNode))
    names(dat) <- c("ALG","node","Count")
    return(dat)
}

all_node_asn <- lapply(paste0(directory_with_nodes, all_node_files), load_node)
names(all_node_asn) <- sapply(strsplit(all_node_files, '_'), function(x){x[1]})
# all_node_asn is a list node -> busco assignments gruped by assigments

######### Processing rearrangement file

all_rearrangements <- read.table(rearrangement_file, sep = '\t')
nodes_with_changes <- list()
nodes_with_changes[names(table(all_rearrangements[, 2]))] <- TRUE

# if(is.null(nodes_with_changes[['n84']])){ print('well hellp')}
# is.null(nodes_with_changes[['n583']])

######### Loading tree
tree <- read.tree(treefile)
# rooted_tree <- tree # root(tree, "Panorpa_germanica", resolve.root = T)
tree_components <- c(tree$tip.label, tree$node.label)

######### Species 2 Family list!

library('gsheet')
library('phytools')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

# species2family <- list()
# species2family[all_genome_data[, 'species']] <- all_genome_data[, 'family']
family2species <- aggregate(all_genome_data[, 'species'], list(all_genome_data[, 'family']), c)
# this looks graet
family2species <- family2species[family2species[, 1] != '', ]
family2node <- list()
family2descendants <- list()

head(family2species)

for (family_i in 1:nrow(family2species)){
    family <- family2species[family_i, 1]
    species <- unlist(family2species[family_i, 2])

    # keep only the species in the tree
    species <- species[species %in% tree$tip.label]

    if (length(species) > 1){
        mrca <- getMRCA(tree, species)
        family2node[[family]] <- tree_components[mrca]
        family2descendants[[family]] <- tree_components[getDescendants(tree, mrca, curr=NULL)]
    }
}

nodes_to_avoid <- unique(unlist(family2descendants))


########## Plotting the tree

pdf(output, height = 60, width = 20)

plot(tree,  show.node.label = TRUE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

chwidth <- max(lastPP$x.lim) / 100
chheight <- 2
chgap <- max(lastPP$x.lim) / 400
y_smallstep <- 0.04

for (node_i in c(length(tree$tip.label) + 2):length(lastPP$xx)){
    node <- tree_components[node_i]    
    if ( is.null(nodes_with_changes[[node]]) ){
        print(paste("Skipping", node, " - no changes."))
        next
    }
    if ( node %in% nodes_to_avoid){
        print(paste("Skipping", node, " - subfamily node."))
        next
    }
    if ( node %in% unlist(family2node) ){
        print(paste("Skipping", node, " - Family node."))
        next
    }
    node_asn <- all_node_asn[[node]]
    
    x <- lastPP$xx[node_i]
    y <- lastPP$yy[node_i]

    # get order right 
    lgs_to_plot <- levels(node_asn[, 'node'])
    mixed <- 0
    all_node_lgs <- unique(node_asn[, 2])
    dominan_asn <- sapply(all_node_lgs, function(x){ chr_subset = node_asn[node_asn[, 2] == x, ]; chr_subset[which.max(chr_subset[, 'Count']), 1]})
    ordered_lgs <- all_node_lgs[order(dominan_asn, decreasing = T)]
    lgs_to_plot <- ordered_lgs
    ##

    yfrom <- y + y_smallstep
    yto <- y + y_smallstep + chheight
    for (chr in 1:length(lgs_to_plot)){
        xfrom <- x - (chr * (chwidth + chgap))
        xto <- x - ((chr - 1) * (chwidth + chgap) + chgap)

        bar_subset <- node_asn[node_asn[, 'node'] == lgs_to_plot[chr], ]
        # making sure the order in the plot will be consisent
        bar_subset <- bar_subset[order(bar_subset[, 'ALG'], decreasing = F), ]
        if (simplified == 0){
            columns_to_plot <- y + y_smallstep + (chheight * c(0, cumsum(bar_subset[, "Count"])) / sum(bar_subset[, "Count"]))
        } else {
            theshold <- 1 / simplified
            pieces_per_ALG <- round((bar_subset[, "Count"] / sum(bar_subset[, "Count"])) / theshold)
            mixed <- (simplified - sum(pieces_per_ALG))
            columns_to_plot <- y + y_smallstep + (chheight * c(0, cumsum(pieces_per_ALG)) / simplified)
        }

        for(i in 1:nrow(bar_subset)){
            if ( columns_to_plot[i] != columns_to_plot[i + 1]){
                rect(xfrom, columns_to_plot[i], 
                    xto, columns_to_plot[i + 1], 
                    col = pal[as.character(bar_subset[i, 'ALG'])], border = NA)
            }
        }
        if (mixed != 0){
            rect(xfrom, columns_to_plot[i + 1], 
                xto, yto, 
                col = 'gray', border = NA)            
        }
    }
}

dev.off()

familes_to_plot <- rev(read.table('tables/families_ordered_by_tree.tsv')[, 1])

algs_per_family <- sapply(family2node, function(x){length(unique(all_node_asn[[x]][, 'node']))})
max_nodes_to_plot <- max(algs_per_family)
barwidth <- 0.45

pdf('figures/family_paints.pdf', height = 60, width = 10)

plot(NA, xlim = c(-1, max_nodes_to_plot + barwidth), ylim = c(0, length(familes_to_plot)), bty = 'n', axes = F, xlab = '', ylab = '')

for (rowi in 1:length(familes_to_plot)){
    ybaseline <- rowi - 1
    family <- familes_to_plot[rowi]

    text(-0.5, ybaseline, family, pos = 3)

    family_node <- family2node[[family]]
    if ( is.null(family_node) ){
        print(paste("Skipping", family, " - not in the syngraph tree"))
        next
    }
    dat <- all_node_asn[[family_node]]
    if ( is.null(dat) ) {
        print(paste("Skipping", family, " - lacking family node reconstruction"))
    }

    nodes_to_plot <- levels(dat[, 'node'])
    dominan_asn <- sapply(nodes_to_plot, function(x){ chr_subset = dat[dat[, 2] == x, ]; chr_subset[which.max(chr_subset[, 'Count']), 1]})
    ordered_lgs <- nodes_to_plot[order(dominan_asn, decreasing = T)]
    nodes_to_plot <- rev(ordered_lgs)

    for (node in 1:length(nodes_to_plot)){
        xpos <- node
        bar_subset <- dat[dat[, 'node'] == nodes_to_plot[node], ]
        # making sure the order in the plot will be consisent
        bar_subset <- bar_subset[order(bar_subset[, 'ALG'], decreasing = F), ]
        columns_to_plot <- c(0, cumsum(bar_subset[, "Count"])) / sum(bar_subset[, "Count"]) * 0.8

        for(i in 1:nrow(bar_subset)){
            rect(xpos - barwidth, ybaseline + columns_to_plot[i], 
                xpos + barwidth, ybaseline + columns_to_plot[i + 1], 
                col = pal[as.character(bar_subset[i, 'ALG'])], border = NA)
        }
    }
}

dev.off()