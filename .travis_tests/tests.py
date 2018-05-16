#!/usr/bin/env python
from copy import copy
from gzip import open as gopen
from os.path import dirname,realpath
from treeswift import read_tree_newick
TREESTR = gopen('%s/test.tre.gz'%dirname(realpath(__file__))).read().decode().strip()

# tests
def test_closest_leaf_to_root(t):
    l,d = t.closest_leaf_to_root()
def test_coalescence_waiting_times(t):
    for l in t.coalescence_waiting_times():
        pass
def test_collapse_short_branches(t):
    t.collapse_short_branches(float('inf'))
def test_copy(t):
    o = copy(t)
def test_diameter(t):
    d = t.diameter()
def test_distances_from_root(t):
    for n,d in t.distances_from_root():
        pass
def test_edge_length_sum(t):
    o = t.edge_length_sum()
    o = t.edge_length_sum(leaves=False)
    o = t.edge_length_sum(internal=False)
    o = t.edge_length_sum(leaves=False, internal=False)
def test_extract_tree_with(t):
    pass # TODO
def test_extract_tree_without(t):
    pass # TODO
def test_furthest_from_root(t):
    n,d = t.furthest_from_root()
def test_label_to_node(t):
    for l,n in t.label_to_node().items():
        pass
def test_mrca(t):
    pass # TODO
def test_newick(t):
    s = t.newick(); s = str(t)
def test_resolve_polytomies(t):
    t.resolve_polytomies() # TODO
def test_scale_edges(t):
    t.scale_edges(1.5)
def test_suppress_unifurcations(t):
    t.suppress_unifurcations() # TODO
def test_traverse_inorder(t):
    for n in t.traverse_inorder():
        pass
def test_traverse_leaves(t):
    for l in t.traverse_leaves():
        pass
def test_traverse_levelorder(t):
    for n in t.traverse_levelorder():
        pass
def test_traverse_postorder(t):
    for n in t.traverse_postorder():
        pass
def test_traverse_preorder(t):
    for n in t.traverse_preorder():
        pass
def test_traverse_rootdistorder(t):
    for n in t.traverse_rootdistorder(ascending=True):
        pass
    for n in t.traverse_rootdistorder(ascending=False):
        pass

# run tests
if __name__ == "__main__":
    tests = [v for k,v in locals().items() if callable(v) and v.__module__ == __name__]
    for test in tests:
        test(read_tree_newick(TREESTR))