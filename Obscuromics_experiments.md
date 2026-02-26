
### Scripts

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

scripts/prep_syngraph_tree.py data/hymenoptera/hymenoptera.supermatrix.phy.treefile data/hymenoptera/excluded_syngraph_busco_tables data/hymenoptera/hymenoptera.pruned.NEW.newick

syngraph build -d data/hymenoptera/syngraph_busco_tables -m -o data/hymenoptera/syngraph_run_pruned/hymenoptera.syngraph_build.pruned

syngraph infer -g data/hymenoptera/syngraph_run_pruned/hymenoptera.syngraph_build.pruned.pickle -t data/hymenoptera/hymenoptera.pruned.NEW.newick -m 50 -r 2 -a quick -s Panorpa_germanica -o data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer -d

# success

SYNGTAB=data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_tabulate.table.tsv
syngraph tabulate -g data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.with_ancestors.pickle -o data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_tabulate
python3 scripts/decompose_table.py --syngtab $SYNGTAB --out-dir data/hymenoptera/syngraph_run_pruned2/node_asn

Rscript scripts/plot_tree_with_ALGs_at_nodes.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv -n data/hymenoptera/syngraph_run_pruned2/node_asn/ -l family -s data/hymenoptera/Hymenoptera_genomes_taxonomy.tsv -o figures/syngraph_tree_of_changes_pruned_2

Rscript scripts/markers_per_node.R -t data/hymenoptera/syngraph_run_pruned2/hymenoptera.mindist.m50.syngraph_infer.newick.txt -nodes_dir data/hymenoptera/syngraph_run_pruned2/node_asn/ -a data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.n4.ALGs.tsv  -o figures/syngraph_tree_pruned2_makers_in_nodes


Rscript scripts/define_ALGs.R -o data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv -n n1 -lgn y -i 'data/hymenoptera/syngraph_run_pruned2/node_asn/'

data/hymenoptera/syngraph_run/ALG_syngraph.hymenoptera.prune2.n1.ALGs.tsv
```