#!/usr/bin/env python
from gzip import open as gopen
from os.path import dirname,realpath
from treeswift import read_tree_newick
TREESTR = gopen('%s/test.tre.gz'%dirname(realpath(__file__))).read().decode().strip()

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
    tests = [v for k,v in locals().items() if callable(v) and v.__module__ == __name__]
    for test in tests:
        test(read_tree_newick(TREESTR))