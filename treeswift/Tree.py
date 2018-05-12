#! /usr/bin/env python
from treeswift.Node import Node
from copy import copy
from gzip import open as gopen
from warnings import warn
try:                # Python 3
    from queue import Queue
except ImportError: # Python 2
    from Queue import Queue
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

    def closest_leaf_to_root(self):
        '''Return the leaf that is closest to the root and the corresponding distance. Edges with no length will be considered to have a length of 0

        Returns:
            tuple: First value is the closest leaf to the root, and second value is the corresponding distance
        '''
        best = (None,float('inf')); d = dict()
        for node in self.traverse_preorder():
            d[node] = {True:0,False:node.edge_length}[node.edge_length is None]
            if not node.is_root():
                d[node] += d[node.parent]
            if node.is_leaf() and d[node] < best[1]:
                best = (node,d[node])
        return best

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

    def extract_tree(self, labels, without, suppress_unifurcations=True):
        '''Helper function for extract_tree_* functions'''
        if not isinstance(labels, set):
            labels = set(labels)
        label_to_leaf = dict(); keep = set()
        for node in self.traverse_leaves():
            label_to_leaf[str(node)] = node
            if (without and str(node) not in labels) or (not without and str(node) in labels):
                keep.add(node)
        for node in list(keep):
            for a in node.traverse_ancestors(include_self=False):
                keep.add(a)
        out = Tree(); out.root.label = self.root.label; out.root.edge_length = self.root.edge_length
        q_old = Queue(); q_old.put(self.root)
        q_new = Queue(); q_new.put(out.root)
        while not q_old.empty():
            n_old = q_old.get(); n_new = q_new.get()
            for c_old in n_old.children:
                if c_old in keep:
                    c_new = Node(label=str(c_old), edge_length=c_old.edge_length); n_new.add_child(c_new)
                    q_old.put(c_old); q_new.put(c_new)
        if suppress_unifurcations:
            out.suppress_unifurcations()
        return out

    def extract_tree_without(self, labels, suppress_unifurcations=True):
        '''Extract a copy of this Tree without the leaves labeled by the strings in `labels`

        Args:
            labels (set): Set of leaf labels to exclude

        Returns:
            Tree: Copy of this Tree, exluding the leaves labeled by the strings in `labels`
        '''
        return self.extract_tree(labels, True, suppress_unifurcations)

    def extract_tree_with(self, labels, suppress_unifurcations=True):
        '''Extract a copy of this Tree with only the leaves labeled by the strings in `labels`

        Args:
            leaves (set): Set of leaf labels to include.

        Returns:
            Tree: Copy of this Tree, including only the leaves labeled by the strings in `labels`
        '''
        return self.extract_tree(labels, False, suppress_unifurcations)

    def furthest_from_root(self):
        '''Return the Node that is furthest from the root and the corresponding distance. Edges with no length will be considered to have a length of 0

        Returns:
            tuple: First value is the furthest Node from the root, and second value is the corresponding distance
        '''
        best = (self.root,0); d = dict()
        for node in self.traverse_preorder():
            d[node] = {True:0,False:node.edge_length}[node.edge_length is None]
            if not node.is_root():
                d[node] += d[node.parent]
            if d[node] > best[1]:
                best = (node,d[node])
        return best

    def get_nodes_with_label(self, labels):
        '''Return a dictionary with all nodes labeled by a label in `labels`. If multiple nodes are labeled by a given label, only the last (preorder traversal) will be obtained

        Args:
            labels (set): Set of leaf labels to get

        Returns:
            dict: Dictionary mapping labels to the corresponding nodes
        '''
        if not isinstance(labels, set):
            labels = set(labels)
        label_to_node = dict()
        for node in self.traverse_preorder():
            if str(node) in labels and str(node):
                label_to_node[str(node)] = node
        if len(label_to_node) != len(labels):
            warn("Not all given labels exist in the tree")
        return label_to_node

    def mrca(self, labels):
        '''Return the Node that is the MRCA of the nodes labeled by a label in `labels`. If multiple nodes are labeled by a given label, only the last (preorder traversal) will be obtained

        Args:
            labels (set): Set of leaf labels

        Returns:
            Node: The MRCA of the Node objects labeled by a label in `labels`
        '''
        label_to_node = self.get_nodes_with_label(labels)
        count = dict()
        for node in label_to_node.values():
            for a in node.traverse_ancestors():
                if a not in count:
                    count[a] = 0
                count[a] += 1
                if count[a] == len(label_to_node):
                    return a
        assert False, "There somehow does not exist an MRCA for the given labels"

    def newick(self):
        '''Output this Tree as a Newick string

        Returns:
            str: Newick string of this Tree
        '''
        if self.root.edge_length is None:
            return '%s;' % self.root.newick()
        else:
            return '%s:%f;' % (self.root.newick(), self.root.edge_length)

    def suppress_unifurcations(self):
        '''Remove all nodes with only one child and directly attach child to parent'''
        q = Queue(); q.put(self.root)
        while not q.empty():
            node = q.get()
            if len(node.children) != 1:
                for c in node.children:
                    q.put(c)
                continue
            child = node.children.pop()
            if node.is_root():
                self.root = child; child.parent = None
            else:
                parent = node.parent; parent.remove_child(node); parent.add_child(child)
            if node.edge_length is not None:
                if child.edge_length is None:
                    child.edge_length = 0
                child.edge_length += node.edge_length
            q.put(child)

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
