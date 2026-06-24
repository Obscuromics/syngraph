suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-i", "-input",  
    dest="i", help="Name of the input directory with all the nodes assignments (names nXXX_seq.tsv).")
parser$add_argument("-o", "-output",
    dest="o", help="Name of the output table with defined ALGs 1-5.")
parser$add_argument("-n", "-node", default = 'n1',  
    dest="n", help="Spedify node(s) that defines the linkage groups (default n1), if two nodes separate by _ e.g. 'n1_n2' to indicate combined assignment from n1 and n2.")
parser$add_argument("-r", "-reference", default = '',
    dest="r", help="Specify reference assignments to define names (not assignments) for the ALGs. (e.g. to keep them consistent with previous definitions).")
parser$add_argument("-as", "-alg-size", default = 0.05,  
    dest="as", help="Specify proportion of assigned BUSCOs for the ALGs to be considered. (relevant for assignments of two nodes)", type = 'double')
parser$add_argument("-lgn", "-linkage_group_name", default = 'd',  
    dest="lgn", help="The prefix name for the specified lg")

args <- parser$parse_args()
# Rscript scripts/define_ALGs.R -o data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv -n n4 -lgn y -i 'data/hymenoptera/syngraph_run/node_asn/'
# args$n <- 'n4'
# args$n <- 'n4_n5'
# args$o <- 'data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv'
# args$lgn <- 'y'
# args.r <- NA
# args$i <- 'data/hymenoptera/syngraph_run/node_asn/'
# args$as <- 0.05

# read bibio BUSCOs;
# read old assignments
if( args$r != '' ){
    ref_algs <- read.table(args$r, col.names = c('busco', 'ALG'))
    row.names(ref_algs) <- ref_algs[, 'busco']
}

syngraph <- read.table(paste0(args$i, "/busco_order.tsv"), col.names = c('busco'))
alg_groups <- list()

if ("_" %in% strsplit(args$n, "")[[1]]){
    n1_name <- strsplit(args$n, "_")[[1]][1]
    n2_name <- strsplit(args$n, "_")[[1]][2]
    # read n1
    # read n2

    syngraph[, n1_name] <- read.table(paste0(args$i, "/", n1_name, "_seq.tsv"), blank.lines.skip = F)[, 1]
    syngraph[, n2_name] <- read.table(paste0(args$i, "/", n2_name, "_seq.tsv"), blank.lines.skip = F)[, 1]
    syngraph[, 'node'] <- paste(syngraph[, n1_name], syngraph[, n2_name])

    # keep only those with assignments in both
    syngraph <- syngraph[syngraph[, n1_name] != '' & syngraph[, n2_name] != '', ]
    syngraph[, 'ALG'] <- NA

    # 
    print(paste("ALGs in ", n1_name, ":", length(unique(syngraph[, n1_name]))))
    print(paste("ALGs in ", n2_name, ":", length(unique(syngraph[, n2_name]))))

    # sorted table shows 5 corresponding ALGs + several genes that moved around, those we can ignore
    table_of_groups <- sort(table(paste(syngraph[, n1_name], syngraph[, n2_name])), T) # 4 LGs, not 5
    rel_sizes_of_groups <- as.numeric(table_of_groups) / sum(as.numeric(table_of_groups))
    number_of_ALGs <- sum(rel_sizes_of_groups > args$as)
    print(table_of_groups)
    data.frame(group = names(table_of_groups), size = as.numeric(table_of_groups), rel_size = as.numeric(rel_sizes_of_groups), ALG_pass = rel_sizes_of_groups > args$as)  

    dominant_alg_pairs <- strsplit(names(table_of_groups)[1:number_of_ALGs], ' ')

    for (ALG in 1:number_of_ALGs){
        n1_lg <- dominant_alg_pairs[[ALG]][1]
        n2_lg <- dominant_alg_pairs[[ALG]][2]
        alg_rows <- which(syngraph[, n1_name] == n1_lg & syngraph[, n2_name] == n2_lg)

        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }
} else {
    syngraph[, 'node'] <- read.table(paste0(args$i, "/", args$n, "_seq.tsv"), blank.lines.skip = F)[, 1]

    syngraph <- syngraph[syngraph[, 'node'] != '', ]
    syngraph[, 'ALG'] <- NA

    table_of_groups <- sort(table(syngraph[, 'node']), T)
    number_of_ALGs <- length(table_of_groups)
    print(table_of_groups)

    for ( ALG in 1:length(table_of_groups)){
        alg_name <- names(table_of_groups)[ALG]
        alg_rows <- which(syngraph[, 'node'] == alg_name)
        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }

}

row.names(syngraph) <- syngraph[, 'busco']
node2lg <- paste0(args$lgn, 1:number_of_ALGs)
names(node2lg) <- names(table_of_groups)[1:number_of_ALGs]

# something like this needs to be done when reference is provided
# for (ALG in 1:length(table_of_groups)){
#     node_name <- names(table_of_groups)[ALG]
#     ALG2Bibio <- sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)
#     chr_in_bibio <- names(ALG2Bibio[1])
#     alg_label <- paste0(args$lgn, bibio_odb10_asn[chr_in_bibio, 'alg'])
#     print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))

#     if ( alg_label %in% node2lg){
#         node2lg[alg_label == node2lg] <- paste0(alg_label, 'a')
#         node2lg[node_name] <- paste0(alg_label, 'b')
#     } else {
#         node2lg[node_name] <- alg_label
#     }
# }

syngraph[, 'ALG'] <- node2lg[syngraph[, 'node']]
    # syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label

final_ALGs <- rbind(syngraph[!is.na(syngraph[, 'ALG']), c('busco', 'ALG')])
final_ALGs <- final_ALGs[order(final_ALGs[, 2], final_ALGs[, 1]), ]

# 'data/syngraph.pruned2.100.ALGs.inferred.tsv'
# args$o = 'tables/ALG_syngraph.pruned2.100.ALGs_brachycera.tsv'
write.table(final_ALGs, file = args$o, col.names = F, row.names = F, quote = F, sep = '\t')