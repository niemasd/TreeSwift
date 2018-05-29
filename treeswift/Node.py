#! /usr/bin/env python
from collections import deque
from copy import copy
try:                # Python 3
    from queue import PriorityQueue
except ImportError: # Python 2
    from Queue import PriorityQueue
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
        self.children = list()         # list of child Node objects
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
        if not isinstance(child, Node):
            raise TypeError("child must be a Node")
        self.children.append(child); child.parent = self

    def child_nodes(self):
        '''Return a `list` containing this Node object's children

        Returns:
            set: A `list` containing this Node object's children
        '''
        return copy(self.children)

    def contract(self):
        '''Contract this `Node` by directly connecting its children to its parent'''
        if self.is_root():
            return
        for c in self.children:
            if self.edge_length is not None and c.edge_length is not None:
                c.edge_length += self.edge_length
            self.parent.add_child(c)
        self.parent.remove_child(self)

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
        if not isinstance(child, Node):
            raise TypeError("child must be a Node")
        try:
            self.children.remove(child); child.parent = None
        except:
            raise RuntimeError("Attempting to remove non-existent child")

    def traverse_ancestors(self, include_self=True):
        '''Traverse over the ancestors of this Node

        Args:
            include_self (bool): True to include self in the traversal, otherwise False
        '''
        if not isinstance(include_self, bool):
            raise TypeError("include_self must be a bool")
        if include_self:
            c = self
        else:
            c = self.parent
        while c is not None:
            yield c; c = c.parent

    def traverse_inorder(self):
        '''Perform an inorder traversal starting at this Node object'''
        c = self; s = deque(); done = False
        while not done:
            if c is None:
                if len(s) == 0:
                    done = True
                else:
                    c = s.pop(); yield c
                    if len(c.children) == 0:
                        c = None
                    elif len(c.children) == 2:
                        c = c.children[1]
                    else:
                        raise RuntimeError(INORDER_NONBINARY)
            else:
                s.append(c)
                if len(c.children) == 0:
                    c = None
                elif len(c.children) == 2:
                    c = c.children[0]
                else:
                    raise RuntimeError(INORDER_NONBINARY)

    def traverse_internal(self):
        '''Traverse over the internal nodes below (and including) this Node object'''
        for n in self.traverse_preorder():
            if not n.is_leaf():
                yield n

    def traverse_leaves(self):
        '''Traverse over the leaves below this Node object'''
        for n in self.traverse_preorder():
            if n.is_leaf():
                yield n

    def traverse_levelorder(self):
        '''Perform a levelorder traversal starting at this Node object'''
        q = deque(); q.append(self)
        while len(q) != 0:
            n = q.popleft(); yield n
            for c in n.children:
                q.append(c)

    def traverse_postorder(self):
        '''Perform a postorder traversal starting at this Node object'''
        s1 = deque(); s2 = deque(); s1.append(self)
        while len(s1) != 0:
            n = s1.pop(); s2.append(n)
            for c in n.children:
                s1.append(c)
        while len(s2) != 0:
            yield s2.pop()

    def traverse_preorder(self):
        '''Perform a preorder traversal starting at this Node object'''
        s = deque(); s.append(self)
        while len(s) != 0:
            n = s.pop(); yield n
            for c in n.children:
                s.append(c)

    def traverse_rootdistorder(self, ascending=True):
        '''Perform a traversal of the Node objects in the subtree rooted at this Node in either ascending (`ascending=True`) or descending (`ascending=False`) order of distance from this Node'''
        if not isinstance(ascending, bool):
            raise TypeError("ascending must be a bool")
        pq = PriorityQueue(); dist_from_root = dict()
        for node in self.traverse_preorder():
            if node == self:
                d = 0
            else:
                d = dist_from_root[node.parent] + node.edge_length
            dist_from_root[node] = d
            if ascending:
                pq.put((d,node))
            else:
                pq.put((-d,node))
        while not pq.empty():
            priority,node = pq.get()
            if ascending:
                yield (priority,node)
            else:
                yield (-priority,node)