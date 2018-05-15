#! /usr/bin/env python
from copy import copy
try:                # Python 3
    from queue import Queue
except ImportError: # Python 2
    from Queue import Queue
INORDER_NONBINARY = "Can't do inorder traversal on non-binary tree"
INVALID_NEWICK = "Tree not valid Newick tree"

class Node:
    '''Node class'''
    def __init__(self, label=None, edge_length=None):
        '''Node constructor

        Args:
            label (str): Label of this Node
            edge_length (float): Length of the edge incident to this Node

        Returns:
            Node object
        '''
        self.children = set()          # set of child Node objects
        self.parent = None             # parent Node object (None for root)
        self.label = label             # label
        self.edge_length = edge_length # length of incident edge

    def __lt__(self, other):
        '''Less Than operator. Just compares labels'''
        if not isinstance(other,Node):
            raise TypeError("'<' not supported between instances of '%s' and '%s'"%(type(self).__name__,type(other).__name__))
        elif self.label is None and other.label is not None:
            return True
        elif other.label is None:
            return False
        return self.label < other.label

    def __str__(self):
        '''Represent Node as a string

        Returns:
            str: string representation of this Node
        '''
        return {True:'',False:self.label}[self.label is None]

    def __copy__(self):
        '''Copy this Node

        Returns:
            Node: A copy of this Node
        '''
        out = Node(label=copy(self.label), edge_length=copy(self))
        out.children = copy(self.children)
        out.parent = self.parent
        return out

    def add_child(self, child):
        '''Add child to Node object

        Args:
            child (Node): The child Node to be added
        '''
        if child in self.children:
            raise RuntimeError("Attempting to add existing child")
        self.children.add(child); child.parent = self

    def child_nodes(self):
        '''Return a `set` containing this Node object's children

        Returns:
            set: A `set` containing this Node object's children
        '''
        return copy(self.children)

    def is_leaf(self):
        '''Returns True if this is a leaf

        Returns:
            bool: True if this is a leaf, otherwise False
        '''
        return len(self.children) == 0

    def is_root(self):
        '''Returns True if this is the root

        Returns:
            bool: True if this is the root, otherwise False
        '''
        return self.parent is None

    def newick(self):
        '''Recursive Newick string conversion starting at this Node object

        Returns:
            str: Recursive Newick string conversion starting at this Node object
        '''
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
        '''Remove child from Node object

        Args:
            child (Node): The child to remove
        '''
        if child not in self.children:
            raise RuntimeError("Attempting to remove non-existent child")
        self.children.remove(child); child.parent = None

    def traverse_ancestors(self, include_self=True):
        '''Traverse over the ancestors of this Node

        Args:
            include_self (bool): True to include self in the traversal, otherwise False
        '''
        curr = {True:self,False:self.parent}[include_self]
        while curr is not None:
            yield curr; curr = curr.parent

    def traverse_inorder(self):
        '''Perform an inorder traversal starting at this Node object'''
        c = list(self.children)
        if len(c) != 0 and len(c) != 2:
            raise RuntimeError(INORDER_NONBINARY)
        if len(c) != 0:
            for y in c[0].traverse_inorder():
                yield y
        yield self
        if len(c) != 0:
            for y in c[1].traverse_inorder():
                yield y

    def traverse_internal(self):
        '''Traverse over the internal nodes below (and including) this Node object'''
        for n in self.traverse_levelorder(leaves=False):
            yield n

    def traverse_leaves(self):
        '''Traverse over the leaves below this Node object'''
        for n in self.traverse_levelorder(internal=False):
            yield n

    def traverse_levelorder(self, leaves=True, internal=True):
        '''Perform a levelorder traversal starting at this Node object'''
        q = Queue(); q.put(self)
        while not q.empty():
            n = q.get()
            if n.is_leaf():
                if leaves:
                    yield n
            elif internal:
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
