#! /usr/bin/env python
from treeswift.Node import Node
from copy import copy
from gzip import open as gopen
INVALID_NEWICK = "Tree not valid Newick tree"

class Tree:
    '''Tree class'''
    def __init__(self):
        '''Tree constructor'''
        self.root = Node()  # root Node object

    def __str__(self):
        '''Represent this Tree as a string

        Returns:
            str: string representation of this Tree (Newick string)
        '''
        return self.newick()

    def diameter(self):
        '''Compute the diameter (maximum leaf pairwise distance) of this Tree

        Returns:
            float: The diameter of this Tree
        '''
        d = dict(); best = float('-inf')
        for node in self.traverse_postorder():
            if node.is_leaf():
                d[node] = 0
            else:
                dists = sorted(d[c]+c.edge_length for c in node.children)
                d[node] = dists[-1]; max_pair = dists[-1]+dists[-2]
                if max_pair > best:
                    best = max_pair
        return best

    def distances_from_root(self, leaves=True, internal=True):
        '''Generator over the root-to-tip distances of this Tree'''
        if leaves or internal:
            d = dict()
            for node in self.traverse_preorder():
                if node.is_root():
                    d[node] = {True:0,False:node.edge_length}[node.edge_length is None]
                else:
                    d[node] = d[node.parent] + node.edge_length
                if (node.is_leaf() and leaves) or (not node.is_leaf() and internal):
                    yield d[node]

    def newick(self):
        '''Output this Tree as a Newick string

        Returns:
            str: Newick string of this Tree
        '''
        if self.root.edge_length is None:
            return '%s;' % self.root.newick()
        else:
            return '%s:%f;' % (self.root.newick(), self.root.edge_length)

    def traverse_inorder(self):
        '''Perform an inorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_inorder():
            yield node

    def traverse_internal(self):
        '''Traverse over the internal nodes of this Tree'''
        for node in self.root.traverse_internal():
            yield node

    def traverse_leaves(self):
        '''Traverse over the leaves of this Tree'''
        for node in self.root.traverse_leaves():
            yield node

    def traverse_levelorder(self):
        '''Perform a levelorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_levelorder():
            yield node

    def traverse_postorder(self):
        '''Perform a postorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_postorder():
            yield node

    def traverse_preorder(self):
        '''Perform a preorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_preorder():
            yield node

def read_tree_newick(tree_string):
    '''Read a tree from a Newick string

    Args:
        tree_string (str): Newick string

    Returns:
        Tree: The tree represented by `tree_string`
    '''
    ts = tree_string; t = Tree(); n = t.root
    i = 0
    while i < len(ts):
        if ts[i] == ';':
            assert i == len(ts)-1 and n == t.root, INVALID_NEWICK
        elif ts[i] == '(':
            c = Node(); n.add_child(c); n = c
        elif ts[i] == ')':
            n = n.parent
        elif ts[i] == ',':
            n = n.parent; c = Node(); n.add_child(c); n = c
        elif ts[i] == ':':
            i += 1; ls = ''
            while ts[i] not in {',', ')', ';'}:
                ls += ts[i]; i += 1
            n.edge_length = float(ls); i -= 1
        else:
            label = ''
            while ts[i] not in {':', ',', ';',')'}:
                label += ts[i]; i += 1
            i -= 1; n.label = label
        i += 1
    return t
