suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-i", "-input",  
    dest="i", help="Name of the input directory with all the nodes assignments (names nXXX_seq.table.tsv).")
parser$add_argument("-o", "-output",
    dest="o", help="Name of the output table with defined ALGs 1-5.")
parser$add_argument("-n", "-node", default = 'n1',  
    dest="n", help="Spedify node(s) that defines the linkage groups (default n1), if two nodes separate by _ e.g. 'n1_n2' to indicate combined assignment from n1 and n2.")
parser$add_argument("-r", "-reference", default = NA,  
    dest="r", help="Specify reference assignments to define names (not assignments) for the ALGs. (e.g. to keep them consistent with previous definitions).")
parser$add_argument("-as", "-alg-size", default = 0.05,  
    dest="as", help="Specify proportion of assigned BUSCOs for the ALGs to be considered. (relevant for assignments of two nodes)", type = 'double')
parser$add_argument("-lgn", "-linkage_group_name", default = 'd',  
    dest="lgn", help="The prefix name for the specified lg")

args <- parser$parse_args()
# Rscript scripts/20250719_ALG_definition.R -o tables/ALG_syngraph.m165.ALGs_schizophora_syrphidae.tsv -n n133 -lgn ds

# args$o <- 'tables/ALG_brachycera'
# args$n <- 'n21'
# args$lgn <- 'db'

# args$o <- 'tables/ALG_syngraph.m165.ALGs_schizophora_syrphidae.tsv'
# args$n <- 'n133'
# args$lgn <- 'ds'

# d - diptra; n1n2
# db - brachycera; n21 (should be n13, but I want all 4 splits)
# ds - schizophora and surphidae; n133
# dm - culicidae (m for mosquitos); n5
# dc - chironomidae -> impossible with the current run;

# read bibio BUSCOs;
# read old assignments
if( !is.na(args.r) ){
    ref_algs <- read.table(args.r, col.names = c('busco', 'ALG'))
    row.names(ref_algs) <- ref_algs[, 'busco']
}

alg_groups <- list()

if ("_" %in% strsplit(args$n, "")[[1]]){
    n1_name <- strsplit(args$n, "_")[[1]][1]
    n2_name <- strsplit(args$n, "_")[[1]][2]
    # read n1
    # read n2
    n1 <- read.table(paste0(args$i, "/", n1_name, "_seq.table.tsv"), col.names = c('busco', n1_name))
    n2 <- read.table(paste0(args$i, "/", n2_name, "_seq.table.tsv"), col.names = c('busco', n2_name))
    syngraph <- merge(n1, n2) # this will drop all those missing in one or the other
    syngraph[, 'ALG'] <- NA

    # double check both those are 5
    print(length(unique(n1[, n1_name])))
    print(length(unique(n2[, n2_name])))

    # sorted table shows 5 corresponding ALGs + several genes that moved around, those we can ignore
    table_of_groups <- sort(table(paste(syngraph[, n1_name], syngraph[, n2_name])), T) # 4 LGs, not 5
    number_of_ALGs <- (table_of_groups / sum(table_of_groups)) > 0.05
    print(table_of_groups)

    dominant_alg_pairs <- strsplit(names(table_of_groups)[1:5], ' ')


    for (ALG in 1:5){
        n1_name <- dominant_alg_pairs[[ALG]][1]
        n2_name <- dominant_alg_pairs[[ALG]][2]
        alg_rows <- which(syngraph[, 'n1'] == n1_name & syngraph[, 'n2'] == n2_name)

        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }
} else {
    syngraph <- read.table(paste0("data/syngraph/node_assignments/", args$n ,"_asgn.tsv"), col.names = c('busco', 'node'))
    syngraph[, 'ALG'] <- NA

    table_of_groups <- sort(table(syngraph[, 'node']), T)
    print(table_of_groups)

    for ( ALG in 1:length(table_of_groups)){
        alg_name <- names(table_of_groups)[ALG]
        alg_rows <- which(syngraph[, 'node'] == alg_name)
        alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    }

}

row.names(syngraph) <- syngraph[, 'busco']
node2lg <- rep('X', length(table_of_groups))
names(node2lg) <- names(table_of_groups)

if (args$lgn == 'd'){
    for (ALG in 1:5){
        chr_in_bibio <- names(sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)[1])
        alg_label <- paste0('d', bibio_odb10_asn[chr_in_bibio, 'alg'])
        
        syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label
        print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))
    }
} 
if (args$lgn == 'ds') {
    for (ALG in 1:length(table_of_groups)){
        node_name <- names(table_of_groups)[ALG]
        ALG2ME <- sort(table(dros_odb12[dros_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)
        chr_in_dros <- names(ALG2ME[1])
        alg_label <- drosophila_chromosomes[chr_in_dros, 'alg']
        print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))
        node2lg[node_name] <- alg_label
    }
} else {
    for (ALG in 1:length(table_of_groups)){
        node_name <- names(table_of_groups)[ALG]
        ALG2Bibio <- sort(table(bibio_odb12[bibio_odb12[, 'busco_odb12'] %in% alg_groups[[ALG]], 'chr']), T)
        chr_in_bibio <- names(ALG2Bibio[1])
        alg_label <- paste0(args$lgn, bibio_odb10_asn[chr_in_bibio, 'alg'])
        print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))

        if ( alg_label %in% node2lg){
            node2lg[alg_label == node2lg] <- paste0(alg_label, 'a')
            node2lg[node_name] <- paste0(alg_label, 'b')
        } else {
            node2lg[node_name] <- alg_label
        }
    }
}

syngraph[, 'ALG'] <- node2lg[syngraph[, 'node']]
    # syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label



### Comment out the next few lines if not to include the dot.
ALGs_pruned_with_dot <- read.table('data/ALG6_BUSCOs.tsv', col.names = c('busco', 'ALG'))
ALG6_buscos <- ALGs_pruned_with_dot[ALGs_pruned_with_dot[, 'ALG'] == 'd6', ] # take only LG from the first good reconstruction of the ancestral state
# sum((ALG6_buscos[, 'busco'] %in% syngraph[, 'busco']))
ALG6_buscos <- ALG6_buscos[!(ALG6_buscos[, 'busco'] %in% syngraph[, 'busco']), ] # remove buscos assigned to other LGs
if ( args$lgn != 'd'){
    ALG6_buscos$ALG <- paste0(args$lgn, 6)
}

final_ALGs <- rbind(syngraph[!is.na(syngraph[, 'ALG']), c('busco', 'ALG')], ALG6_buscos)
final_ALGs <- final_ALGs[order(final_ALGs[, 2], final_ALGs[, 1]), ]

# 'data/syngraph.pruned2.100.ALGs.inferred.tsv'
# args$o = 'tables/ALG_syngraph.pruned2.100.ALGs_brachycera.tsv'
write.table(final_ALGs, file = args$o, col.names = F, row.names = F, quote = F, sep = '\t')