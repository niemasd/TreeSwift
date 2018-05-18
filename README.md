# TreeSwift&nbsp;&nbsp;&nbsp;[![Build Status](https://travis-ci.org/niemasd/TreeSwift.svg?branch=master)](https://travis-ci.org/niemasd/TreeSwift)
TreeSwift is a Python library for parsing, manipulating, and iterating over (rooted) tree structures. TreeSwift places an emphasis on speed.

## Installation
TreeSwift can be installed using `pip`:

```bash
sudo pip install treeswift
```

If you are using a machine on which you lack administrative powers, TreeSwift can be installed locally using `pip`:

```bash
pip install --user treeswift
```

## Usage
Typical usage should be as follows:

1. Import the `treeswift` package
2. Use `treeswift.read_tree_newick` to load your Newick tree
3. Use the various `Tree` class functions on the resulting object as you need

```python
import treeswift
tree = treeswift.read_tree_newick(my_newick_string)
for node in tree.traverse_postorder():
    print(node)
```

Full documentation can be found at [https://niemasd.github.io/TreeSwift/](https://niemasd.github.io/TreeSwift/), and more examples can be found in the [TreeSwift Wiki](https://github.com/niemasd/TreeSwift/wiki).
