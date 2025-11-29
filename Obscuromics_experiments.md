
### Scripts

Adding these scripts

```
scripts/define_ALGs.R
scripts/plot_dotplot.R
scripts/plot_node.R
scripts/plot_syngraph_tree_with_nodes.R
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