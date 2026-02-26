### Requirements

Create conda environment for syngraph:
```
conda create -n syngraph python=3.9 && conda activate syngraph
```

Download obscuromics group modified version of syngraph (Mackintosh et al., 2023):

```
git clone -b desire-path https://github.com/Obscuromics/syngraph.git
cd syngraph
pip install .
```

if pygraphviz fails to build on MacOS because graphviz is missing you can install graphviz with conda (or sudo or brew) and retry:

```
conda install graphviz
pip install .
```

If you would like to run syngraph in the notebook install jupyter to your syngraph conda environment:

```
conda install jupyter
jupyter notebook
```
