#! /usr/bin/env python
'''TreeSwift: Fast tree module for Python 2 and 3'''
INORDER_NONBINARY = "Can't do inorder traversal on non-binary tree"
INVALID_NEWICK = "Tree not valid Newick tree"
from copy import copy
from gzip import open as gopen
try:                # Python 3
    from queue import Queue
except ImportError: # Python 2
    from Queue import Queue
def valid_order(o):
    return o == 'level' or o == 'in' or o == 'post' or o == 'pre'

class Node:
    '''Node class'''
    def __init__(self, parent=None, label=None, edge_length=None, edge_label=None):
        '''Node constructor'''
        self.children = set()          # set of child Node objects
        self.parent = parent           # parent Node object (None for root)
        self.label = label             # label
        self.edge_length = edge_length # length of incident edge

    def __str__(self):
        '''Represent Node as string'''
        return self.label

    def add_child(self, child):
        '''Add child to Node object'''
        assert child not in self.children, "Attempting to add existing child"
        self.children.add(child); child.parent = self

    def child_nodes(self):
        '''Return a set containing this Node object's child Node objects'''
        return copy(self.children)

    def is_leaf(self):
        '''Returns True if this is a leaf'''
        return len(self.children) == 0

    def is_root(self):
        '''Returns True if this is the root'''
        return self.parent is None

    def newick(self):
        '''Recursive Newick string conversion starting at this Node object'''
        if self.is_leaf():
            if self.label is None:
                return ''
            else:
                return self.label
        else:
            out = ['(']
            for c in self.children:
                out.append(c.newick())
                if c.edge_length is not None:
                    out.append(':%f' % c.edge_length)
                out.append(',')
            out.pop() # trailing comma
            out.append(')')
            if self.label is not None:
                out.append(self.label)
            return ''.join(out)

    def remove_child(self, child):
        '''Remove child from Node object'''
        assert child in self.children, "Attempting to remove non-existent child"
        self.children.remove(child); child.parent = None

    def traverse_inorder(self):
        '''Perform an inorder traversal starting at this Node object'''
        c = list(self.children)
        assert len(c) in {0,2}, INORDER_NONBINARY
        if len(c) != 0:
            for y in c[0].traverse_inorder():
                yield y
        yield self
        if len(c) != 0:
            for y in c[1].traverse_inorder():
                yield y

    def traverse_leaves(self):
        '''Traverse over the leaves below this Node object'''
        for n in self.traverse_levelorder(only_leaves=True):
            yield n

    def traverse_levelorder(self, only_leaves=False):
        '''Perform a levelorder traversal starting at this Node object'''
        q = Queue(); q.put(self)
        while not q.empty():
            n = q.get()
            if n.is_leaf() or not only_leaves:
                yield n
            for c in n.children:
                q.put(c)

    def traverse_postorder(self):
        '''Perform a postorder traversal starting at this Node object'''
        for c in self.children:
            for n in c.traverse_postorder():
                yield n
        yield self

    def traverse_preorder(self):
        '''Perform a preorder traversal starting at this Node object'''
        yield self
        for c in self.children:
            for n in c.traverse_preorder():
                yield n

class Tree:
    '''Tree class'''
    def __init__(self):
        '''Tree constructor'''
        self.root = Node()  # root Node object

    def newick(self):
        if self.root.edge_length is None:
            return '%s;' % self.root.newick()
        else:
            return '%s:%f;' % (self.root.newick(), self.root.edge_length)

    def traverse_inorder(self):
        '''Perform an inorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_inorder():
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
        '''Perform an postorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_postorder():
            yield node

    def traverse_preorder(self):
        '''Perform an preorder traversal of the Node objects in this Tree'''
        for node in self.root.traverse_preorder():
            yield node

def read_tree_newick(tree_string):
    '''Read a tree from a Newick string'''
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