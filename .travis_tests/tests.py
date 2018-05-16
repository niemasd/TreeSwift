#!/usr/bin/env python
from copy import copy
from gzip import open as gopen
from os.path import dirname,realpath
from random import sample
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
def test_distance_matrix(t):
    m = t.distance_matrix()
def test_edge_length_sum(t):
    o = t.edge_length_sum()
    o = t.edge_length_sum(leaves=False)
    o = t.edge_length_sum(internal=False)
    o = t.edge_length_sum(leaves=False, internal=False)
def test_extract_tree_with(t):
    o = t.extract_tree_with(sample([str(l) for l in t.traverse_leaves()],10))
def test_extract_tree_without(t):
    o = t.extract_tree_without(sample([str(l) for l in t.traverse_leaves()],10))
def test_furthest_from_root(t):
    n,d = t.furthest_from_root()
def test_label_to_node(t):
    for l,n in t.label_to_node().items():
        pass
def test_mrca(t):
    o = t.mrca(sample({str(l) for l in t.traverse_leaves()},10))
def test_newick(t):
    s = t.newick(); s = str(t)
def test_num_lineages_at(t):
    o = t.num_lineages_at(0)
    o = t.num_lineages_at(1)
    o = t.num_lineages_at(float('inf'))
def test_resolve_polytomies(t):
    t.collapse_short_branches(float('inf'))
    t.resolve_polytomies()
def test_sackin(t):
    o = t.sackin()
    o = t.sackin(None)
    o = t.sackin('yule')
    o = t.sackin('pda')
def test_scale_edges(t):
    t.scale_edges(1.5)
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
def test_treeness(t):
    o = t.treeness()

# run tests
if __name__ == "__main__":
    tests = [v for k,v in locals().items() if callable(v) and v.__module__ == __name__]
    for test in tests:
        test(read_tree_newick(TREESTR))