# TreeSwift
TreeSwift: Fast tree module for Python 2 and 3

## Installation
You should be able to install TreeSwift using `pip`, e.g.:

```bash
pip install treeswift
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

Full documentation can be found at [http://treeswift.readthedocs.io/en/latest/](http://treeswift.readthedocs.io/en/latest/).