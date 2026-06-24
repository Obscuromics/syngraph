
### Scripts

Busco table conversion by Sam

```bash
process prepare_busco_tables{
        publishDir "${params.outdir}/syngraph_busco_tables/", mode:'copy'

        input:
        tuple val(meta), path(busco_dir), path(chromosome_file)

        output:
        tuple val(meta), path("${meta}.syngraph.buscos.tsv")

        script:

        """
        cut -f 1,3,4,5 ${busco_dir}/run_*/full_table.tsv | sed '/^#/d' | awk '\$2 != ""' > ${meta}.busco.reformatted.tsv
        grep -f ${chromosome_file} ${meta}.busco.reformatted.tsv > ${meta}.syngraph.buscos.tsv 
        remove_blank_lines.sh ${meta}.syngraph.buscos.tsv
        """
}
```

#### already working

Plot tree with syngraph nodes. Example

```bash
Rscript scripts/plot_syngraph_tree_with_nodes.R -i data/hymenoptera/syngraph_run/hymenoptera.mindist.m50.newick.txt -o data/hymenoptera/syngraph_run/hymenoptera.mindist.m50.newick.pdf
```

This plot is useful when trying to orient in nodes. Syngraph labels them in order of traversing the tree, so while the numbers are somewhat indicative, one needs to see the tree of know which node corresponds to what clade.

---

Decomposes tabulated file from syngraph into assigments of individual nodes. I made them one column at times including missing values and saved names of busco markers corresponding to the lines in individual file `busco_order.tsv`, other files are saved as `nXXX_seq.tsv`. Example usage

```bash
SYNGTAB=data/hymenoptera/syngraph_run/hymenoptera.syngraph_tabulate.table.tsv
python3 scripts/decompose_table.py --syngtab $SYNGTAB --out-dir data/hymenoptera/syngraph_run/node_asn2
```

This is meant to be a fast way of replacing this oneliner

```bash
node='n1'
awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $SYNGTAB ` > n1_node.tsv
```

---

Define ALGs by specification of the node(s) to define the ALGs (one needs to specify path to the directory with all node definitions). User is expected to provide prefix for the ALGs to avoid confusions between taxonomic groups. We recommend using letters from https://id.tol.sanger.ac.uk/

Example usage:

```bash
Rscript scripts/define_ALGs.R -o data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv -n n4 -lgn y -i 'data/hymenoptera/syngraph_run/node_asn/'
```

---

Plot tree with rearrangements on top of it. Can be used to plot states at nodes of certain taxonomic depth (e.g. family level and avove) and can plot both all states or only rearrnged states (if rerrangement file `-r` is specified). Requires tree, ALG definitions and states in all nodes. The nodes where inference failed are painted in red.

```bash
Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t data/hymenoptera/syngraph_run/hymenoptera.mindist.m50.newick.txt -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv -n data/hymenoptera/syngraph_run/node_asn/ -l family -s data/hymenoptera/Hymenoptera_genomes_taxonomy.tsv -o figures/syngraph_tree_of_changes
```


---

#### in progress

```bash
Rscript scripts/plot_dotplot.R -d ./data/hymenoptera/syngraph_busco_tables/ -s1 Abia_candens -s2 Arge_ochropus -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv -o data/hymenoptera/syngraph_run/dotplot_example
```

#### not adjusted
Adding these scripts

```
scripts/plot_dotplot.R
scripts/plot_node.R
scripts/plot_tree_with_ALGs_at_nodes.R
scripts/prep_syngraph_tree.py
```

They were built for diptera, so they will be disfunctional for now, but I will fix them one hby one now

### Diptera

```bash
mkdir -p data/results/infer_alg_diptera/syngraph_busco_tables/
mkdir -p data/results/infer_alg_diptera/excluded_syngraph_busco_tables/

scp farm:/data/tol/teams/jaron/lustre/users/se13/diptera_alg/data/results/infer_alg_diptera/syngraph_busco_tables/* data/results/infer_alg_diptera/syngraph_busco_tables/
```

### Hymenoptera

```bash
while read sp; do mv data/hymenoptera/syngraph_busco_tables/$sp* data/hymenoptera/excluded_syngraph_busco_tables/; done < data/hymenoptera/exclude

Rscript scripts/root_tree.R -i data/hymenoptera/hymenoptera.supermatrix.phy.treefile -o data/hymenoptera/hymenoptera.supermatrix.phy_rooted.treefile -r "Panorpa_germanica"

scripts/prep_syngraph_tree.py -intree data/hymenoptera/hymenoptera.supermatrix.phy_rooted.treefile -incl_dir data/hymenoptera/syngraph_busco_tables -outtree data/hymenoptera/hymenoptera.pruned.NEW.newick
#  -root "Panorpa_germanica"
python3 scripts/root_tree_Sam.py data/hymenoptera/hymenoptera.pruned.NEW.newick outgroup data/hymenoptera/hymenoptera.pruned.rooted.NEW.newick

syngraph build -d data/hymenoptera/syngraph_busco_tables -m -o data/hymenoptera/syngraph_run_pruned/hymenoptera.syngraph_build.pruned

syngraph infer -g data/hymenoptera/syngraph_run_pruned/hymenoptera.syngraph_build.pruned.pickle -t data/hymenoptera/hymenoptera.pruned.rooted.NEW.newick -m 50 -r 2 -a quick -s Panorpa_germanica -o data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer -d

# success

SYNGTAB=data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_tabulate.table.tsv
syngraph tabulate -g data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.with_ancestors.pickle -o data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_tabulate
python3 scripts/decompose_table.py --syngtab $SYNGTAB --out-dir data/hymenoptera/syngraph_run_pruned2/node_asn

Rscript scripts/markers_per_node.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -nodes_dir data/hymenoptera/syngraph_run_pruned2/node_asn/ -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv  -o figures/syngraph_tree_pruned2_makers_in_nodes

Rscript scripts/define_ALGs.R -o data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv -n n1 -lgn y -i 'data/hymenoptera/syngraph_run_pruned2/node_asn/'
Rscript scripts/define_ALGs.R -o data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n0.ALGs.tsv -n n0 -lgn y -i 'data/hymenoptera/syngraph_run_pruned2/node_asn/'

Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv -n data/hymenoptera/syngraph_run_pruned2/node_asn/ -l family -s data/hymenoptera/Hymenoptera_genomes_taxonomy.tsv -o figures/syngraph_tree_of_changes_pruned_3

# python3 scripts/decompose_table.py --syngtab ../../diptera-ALGs/data/syngraph/diptera.no_plecia.mindist.m165.with_ancestors_tabulate.table.tsv --out-dir ../../diptera-ALGs/data/syngraph/node_asn
# Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t ../../diptera-ALGs/data/syngraph/diptera.no_plecia.mindist.m165.newick.txt -a ../../diptera-ALGs/tables/ALGs_syngraph_diptera.tsv -n ../../diptera-ALGs/data/syngraph/node_asn/ -l family -s ../../diptera-ALGs/tables/diptera_master_sheet.tsv -o figures/diptera_syngraph_tree_of_changes


```





### Sponges

They use GitHub for freeze of genomes (Michael Paulini's repo) and Gert's repo on the phylogenetic analysis.
The final final phylogeny is the one Gert sent in his email, I think it's this: 2025-12_10_C30_ASG_BUSCOS_supermatrix_trimmed_trimal-gappyout_Modb12.fasta.contree.
They are using braker3 Michael Paulini's annotation for this paper.
Gert is going to giving an update. He has odb12 busco runs, he will share the long tables with me.
David Ngugi - assembling/analysing metagenomes of all the sponges;

TODO in regard to input data:
- get BUSCOs of the outgroups!
- get BUSCOs of the other (few) assemblies

```bash
DIR=data/sponges
ORIG_TREE=data/sponges/2025-12_10_C30_ASG_BUSCOS_supermatrix_trimmed_trimal-gappyout_Modb12.fasta.contree
TREE=$DIR/tree.pruned.newick

cat species_table_gca.tsv | cut -f 2 | cut -f 1 -d '.' > tolids_short
paste species_table_gca.tsv tolids_short | cut -f 1,3,4 | tr ' ' '_' | awk '{print $3"_"$1"\t"$2}' > renaming_file

# cp syngraph_busco_tables_GCA/GCA_965654285.1.syngraph.buscos.tsv syngraph_busco_tables/odBieUhxx1_Biemna_sp._UH-2024.syngraph.buscos.tsv 
while read line; do
    sp=$(echo $line | awk '{print $1}')
    gca=$(echo $line | awk '{print $2}')
    cp syngraph_busco_tables_GCA/$gca.syngraph.buscos.tsv syngraph_busco_tables/$sp.syngraph.buscos.tsv
done < renaming_file
# just the two missing are missing...

scripts/prep_syngraph_tree.py $ORIG_TREE $DIR/syngraph_busco_tables $TREE
mkdir -p data/sponges/syngraph_run/

syngraph build -d $DIR/syngraph_busco_tables -m -o $DIR/syngraph_run/syngraph_build
syngraph infer -g $DIR/syngraph_run/syngraph_build.pickle -t $TREE -m 3 -r 2 -a quick -s ooCorCand1_Corticium_candelabrum -o $DIR/sponge.syngraph_infer -d

syngraph tabulate -g $DIR/sponge.syngraph_infer.with_ancestors.pickle -o $DIR/sponge.syngraph_infer.syngraph_tabulate
python3 scripts/decompose_table.py --syngtab $DIR/sponge.syngraph_infer.syngraph_tabulate.table.tsv --out-dir $DIR/syngraph_run/node_asn

Rscript scripts/define_ALGs.R -o $DIR/syngraph_run/draft_ALGs.tsv -n n1 -lgn p -i $DIR/syngraph_run/node_asn/


Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t $DIR/sponge.syngraph_infer.newick.txt -a $DIR/syngraph_run/draft_ALGs.tsv -n $DIR/syngraph_run/node_asn/ -o figures/sponges_syngraph_tree_of_changes -r $DIR/sponge.syngraph_infer.rearrangements.tsv

Rscript scripts/markers_per_node.R -t $DIR/sponge.syngraph_infer.newick.txt -nodes_dir $DIR/syngraph_run/node_asn/ -a $DIR/syngraph_run/draft_ALGs.tsv  -o figures/sponged_syngraph_makers_in_nodes
```

TODO - mask glass sponges; Many of the groups will be way to far for expected meaningful. One of the large groups

I am starting over again with Gert's BUSCO runs.

```bash
while read line; do
    sp=$(echo $line | awk '{print $1}')
    gca=$(echo $line | awk '{print $2}')
    tolid=$(echo $sp | cut -f 1 -d "_")
    # echo $tolid
    file=sponge_busco_full_tables/"$tolid"*_full_table.tsv

    ls $file
    cut -f 1,3,4,5 $file | sed '/^#/d' | awk '$2 != ""' | grep "SUPER" > syngraph_busco_tables/$sp.syngraph.buscos.tsv 
    # could be based on SUPER
    # rm temp.busco.reformatted.tsv
done < renaming_file

cd ../..
syngraph build -d $DIR/syngraph_busco_tables -m -o $DIR/syngraph_run/syngraph_build
scripts/prep_syngraph_tree.py $ORIG_TREE $DIR/syngraph_busco_tables $TREE

syngraph infer -g $DIR/syngraph_run/syngraph_build.pickle -t $TREE -m 3 -r 2 -a quick -s ooCorCand1_Corticium_candelabrum -o $DIR/sponge.syngraph_infer -d > $DIR/syngraph.log

Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t $DIR/sponge.syngraph_infer.newick.txt -a $DIR/syngraph_run/draft_ALGs.tsv -n $DIR/syngraph_run/node_asn/ -o figures/sponges_syngraph_tree_of_changes -r $DIR/sponge.syngraph_infer.rearrangements.tsv

mkdir -p $DIR/sponges/excluded_busco_tables

mv $DIR/syngraph_busco_tables/ohAphBeat* $DIR/excluded_busco_tables
mv $DIR/syngraph_busco_tables/ohAscMega* $DIR/excluded_busco_tables
mv $DIR/syngraph_busco_tables/ohCorPlag* $DIR/excluded_busco_tables
mv $DIR/syngraph_busco_tables/ohBolCyan* $DIR/excluded_busco_tables

Rscript scripts/markers_per_node.R -t $DIR/sponge.syngraph_infer.newick.txt -nodes_dir $DIR/syngraph_run/node_asn/ -a $DIR/syngraph_run/draft_ALGs.tsv  -o figures/sponged_syngraph_makers_in_nodes

# [+] Inferring median genome for n34 using data from n51, odTheVald1_Thenea_valdiviae, and odCraInfr1_Craniella_infrequens ...
# [=] Generated 21 LMSs containing 276 markers
# [=] A total of 386 markers are not assigned to an LMS
mv $DIR/syngraph_busco_tables/odCraInfr* $DIR/excluded_busco_tables
mv $DIR/syngraph_busco_tables/ocSycCili* $DIR/excluded_busco_tables

# [+] Inferring median genome for n27 using data from n41, odIanBast1_Ianthella_basta, and n43 ...
# [=] Generated 4 LMSs containing 74 markers

# n9            8                                    n17         odChoReni1_Chondrosia_reniformis      62
# n1           48                                     n3                                       n4      63
# n4           13                                     n9                                      n10      64
# n3           35                                     n7                                       n8      67
# n28           3                                    n43            odHalCaeu1_Halisarca_caerulea      70
# n10           5                                    n19                 odDysAvar1_Dysidea_avara      71
# n17           7                                    n27                                      n28      77
# n43           2    odChoAust1_Chondrilla_australiensis         odChoCari5_Chondrilla_caribensis      78
# n2            3                                     n5                ocSycCili1_Sycon_ciliatum      90
# n27           4                                    n41               odIanBast1_Ianthella_basta      94


for sp in odSpoOffi odDysAva odCarFoli odRhoOdor odHalCaeu odChoCari odChoAust ooOscLobu; do
    mv $DIR/syngraph_busco_tables/"$sp"* $DIR/excluded_busco_tables
done

syngraph build -d $DIR/syngraph_busco_tables -m -o $DIR/syngraph_run/syngraph_build
scripts/prep_syngraph_tree.py $ORIG_TREE $DIR/syngraph_busco_tables $TREE

syngraph infer -g $DIR/syngraph_run/syngraph_build.pickle -t $TREE -m 3 -r 2 -a quick -s ooCorCand1_Corticium_candelabrum -o $DIR/sponge.syngraph_infer -d > $DIR/syngraph.log

syngraph tabulate -g $DIR/sponge.syngraph_infer.with_ancestors.pickle -o $DIR/sponge.syngraph_infer.syngraph_tabulate
rm $DIR/syngraph_run/node_asn/*
python3 scripts/decompose_table.py --syngtab $DIR/sponge.syngraph_infer.syngraph_tabulate.table.tsv --out-dir $DIR/syngraph_run/node_asn

Rscript scripts/markers_per_node.R -t $DIR/sponge.syngraph_infer.newick.txt -nodes_dir $DIR/syngraph_run/node_asn/  -o figures/sponged_syngraph_makers_in_nodes
Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t $DIR/sponge.syngraph_infer.newick.txt -a $DIR/syngraph_run/draft_ALGs.tsv -n $DIR/syngraph_run/node_asn/ -o figures/sponges_syngraph_tree_of_changes -r $DIR/sponge.syngraph_infer.rearrangements.tsv

Rscript scripts/define_ALGs.R -o $DIR/syngraph_run/draft_ALGs.tsv -n n1 -lgn p -i $DIR/syngraph_run/node_asn/

Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t $DIR/sponge.syngraph_infer.newick.txt -a $DIR/syngraph_run/draft_ALGs.tsv -n $DIR/syngraph_run/node_asn/ -o figures/sponges_syngraph_tree_of_changes -r $DIR/sponge.syngraph_infer.rearrangements.tsv

```


#### Tree plotting tests

```bash
Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv -n data/hymenoptera/syngraph_run_pruned2/node_asn/ -l family -s data/hymenoptera/Hymenoptera_genomes_taxonomy.tsv -o figures/syngraph_TEST
```

or to include also paints of species

```bash
Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv -n data/hymenoptera/syngraph_run_pruned2/node_asn/ -l family -s data/hymenoptera/Hymenoptera_genomes_taxonomy.tsv -o figures/syngraph_TEST -species_marker_dir data/hymenoptera/syngraph_busco_tables
```