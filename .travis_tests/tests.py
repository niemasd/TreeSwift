#!/usr/bin/env python
from gzip import open as gopen
from treeswift import read_tree_newick
TREESTR = gopen('test.tre.gz').read().decode().strip()

# tests
def test_inorder(t):
    for n in t.traverse_inorder():
        pass
def test_levelorder(t):
    for n in t.traverse_levelorder():
        pass
def test_postorder(t):
    for n in t.traverse_postorder():
        pass
def test_preorder(t):
    for n in t.traverse_preorder():
        pass

# run tests
if __name__ == "__main__":
    for k,v in locals().items():
        if callable(v) and v.__module__ == __name__:
            v(read_tree_newick(TREESTR))