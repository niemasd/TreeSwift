#!/usr/bin/env python
from copy import copy
from os.path import dirname, realpath
from random import sample
from treeswift import read_tree_newick, read_tree_nexml, read_tree_nexus
PATH = dirname(realpath(__file__))
NEWICK_FILE = '%s/test.tre' % PATH
NEWICK_GZIP = '%s/test.tre.gz' % PATH
NEWICK_STR = open(NEWICK_FILE).read().strip()
NEXUS_FILE = '%s/test.nex' % PATH
NEXUS_GZIP = '%s/test.nex.gz' % PATH
NEXUS_STR = open(NEXUS_FILE).read().strip()
NEXML_FILE = '%s/test.nexml' % PATH
NEXML_GZIP = '%s/test.nexml.gz' % PATH
NEXML_STR = open(NEXML_FILE).read().strip()

# tests
def test_avg_branch_length(t):
    o = t.avg_branch_length()
    o = t.avg_branch_length(terminal=False)
    o = t.avg_branch_length(internal=False)
def test_branch_lengths(t):
    o = t.branch_lengths()
    o = t.branch_lengths(terminal=False)
    o = t.branch_lengths(internal=False)
def test_closest_leaf_to_root(t):
    l,d = t.closest_leaf_to_root()
def test_coalescence_times(t):
    for d in t.coalescence_times():
        pass
    for d in t.coalescence_times(backward=False):
        pass
def test_coalescence_waiting_times(t):
    for l in t.coalescence_waiting_times():
        pass
    for l in t.coalescence_waiting_times(backward=False):
        pass
def test_collapse_short_branches(t):
    copy(t).collapse_short_branches(float('inf'))
def test_contract_low_support(t):
    copy(t).contract_low_support(float('inf'))
def test_copy(t):
    o = copy(t)
def test_diameter(t):
    d = t.diameter()
def test_distance_between(t):
    u,v = list(t.traverse_leaves())[:2]
    d = t.distance_between(u,v)
def test_distance_matrix(t):
    m = t.distance_matrix()
def test_distances_from_parent(t):
    for n,d in t.distances_from_parent():
        pass
    for n,d in t.distances_from_parent(leaves=False):
        pass
    for n,d in t.distances_from_parent(internal=False):
        pass
    for n,d in t.distances_from_parent(unlabeled=True):
        pass
    for n,d in t.distances_from_parent(leaves=False, internal=False, unlabeled=True):
        pass
def test_distances_from_root(t):
    for n,d in t.distances_from_root():
        pass
    for n,d in t.distances_from_root(leaves=False):
        pass
    for n,d in t.distances_from_root(internal=False):
        pass
    for n,d in t.distances_from_root(unlabeled=True):
        pass
    for n,d in t.distances_from_root(leaves=False, internal=False, unlabeled=True):
        pass
def test_edge_length_sum(t):
    o = t.edge_length_sum()
    o = t.edge_length_sum(terminal=False)
    o = t.edge_length_sum(internal=False)
    o = t.edge_length_sum(terminal=False, internal=False)
def test_extract_tree_with(t):
    o = t.extract_tree_with(sample([str(l) for l in t.traverse_leaves()],10))
def test_extract_tree_without(t):
    o = t.extract_tree_without(sample([str(l) for l in t.traverse_leaves()],10))
def test_furthest_from_root(t):
    n,d = t.furthest_from_root()
def test_gamma_statistic(t):
    g = t.gamma_statistic()
def test_get_edge_length(t):
    for n in t.traverse_preorder():
        l = n.get_edge_length()
def test_get_label(t):
    for n in t.traverse_preorder():
        l = n.get_label()
def test_height(t):
    h = t.height()
def test_indent(t):
    s = t.indent()
def test_label_to_node(t):
    for l,n in t.label_to_node().items():
        pass
def test_labels(t):
    for l in t.labels():
        pass
    for l in t.labels(leaves=False):
        pass
    for l in t.labels(internal=False):
        pass
    for l in t.labels(leaves=False,internal=False):
        pass
def test_ladderize(t):
    t.ladderize()
    t.ladderize(ascending=False)
def test_mrca(t):
    o = t.mrca(sample([str(l) for l in t.traverse_leaves()],10))
def test_mrca_matrix(t):
    m = t.mrca_matrix()
def test_newick(t):
    s = t.newick(); s = str(t)
def test_num_lineages_at(t):
    o = t.num_lineages_at(0)
    o = t.num_lineages_at(1)
    o = t.num_lineages_at(float('inf'))
def test_num_nodes(t):
    o = t.num_nodes()
    o = t.num_nodes(leaves=False)
    o = t.num_nodes(internal=False)
    o = t.num_nodes(leaves=False,internal=False)
def test_order(t):
    t.order('edge_length')
    t.order('edge_length_then_label')
    t.order('edge_length_then_label_then_num_descendants')
    t.order('edge_length_then_num_descendants')
    t.order('edge_length_then_num_descendants_then_label')
    t.order('label')
    t.order('label_then_edge_length')
    t.order('label_then_edge_length_then_num_descendants')
    t.order('label_then_num_descendants')
    t.order('label_then_num_descendants_then_edge_length')
    t.order('num_descendants')
    t.order('num_descendants_then_label')
    t.order('num_descendants_then_label_then_edge_length')
    t.order('num_descendants_then_edge_length')
    t.order('num_descendants_then_edge_length_then_label')
def test_rename_nodes_condense(t):
    m = dict()
    for l in t.traverse_leaves():
        m[str(l)] = 'NIEMA'
    t2 = copy(t)
    t2.rename_nodes(m)
    t2.condense()
def test_resolve_polytomies(t):
    t2 = copy(t)
    t2.collapse_short_branches(float('inf'))
    t2.resolve_polytomies()
def test_sackin(t):
    o = t.sackin()
    o = t.sackin(None)
    o = t.sackin('yule')
    o = t.sackin('pda')
def test_scale_edges(t):
    copy(t).scale_edges(1.5)
def test_set_edge_length(t):
    for n in copy(t).traverse_preorder():
        n.set_edge_length(0)
def test_set_label(t):
    for n in copy(t).traverse_preorder():
        n.set_label('NIEMA')
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
def test_write_tree_newick(t):
    t.write_tree_newick('test_write_tree_newick.tre')
    read_tree_newick('test_write_tree_newick.tre')
    t.write_tree_newick('test_write_tree_newick.tre.gz')
    read_tree_newick('test_write_tree_newick.tre.gz')

# run tests
if __name__ == "__main__":
    tests = [v for k,v in locals().items() if callable(v) and v.__module__ == __name__]
    trees = [read_tree_newick(NEWICK_FILE), read_tree_newick(NEWICK_GZIP), read_tree_newick(NEWICK_STR)]
    for nex in [NEXUS_FILE, NEXUS_GZIP, NEXUS_STR]:
        trees += read_tree_nexus(nex).values()
    for t in trees:
        for test in tests:
            test(t)
    for nexml in [NEXML_FILE, NEXML_GZIP, NEXML_STR]:
        read_tree_nexml(nexml)
