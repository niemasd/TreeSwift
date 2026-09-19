#! /usr/bin/env python
from collections import deque
from copy import copy
from io import StringIO
UNSAFE_SYMBOLS = {';', '(', ')', ',', '[', ']', ':', "'"}
INORDER_NONBINARY = "Can't do inorder traversal on non-binary tree"
INVALID_NEWICK = "Tree not valid Newick tree"

def params_str(params):
    '''Return a string representation of node or edge parameters

    Args:
        ``params`` (``dict``): The node or edge parameters

    Returns:
        ``str``: A string representation of ``params``
    '''
    if isinstance(params, dict):
        return '[%s]' % ','.join(f'{k}={v}' for k, v in params.items())
    else:
        if not isinstance(params, str):
            params = str(params)
        if not (params.startswith('[') and params.endswith(']')):
            params = f'[{params}]'
        return params

def parse_params_str(params_str):
    '''Parse a string representation of node or edge parameters

    Args:
        ``params_str`` (``str``): The string representation of node or edge parameters

    Returns:
        ``dict``: The parsed version of ``params_str``
    '''
    if isinstance(params_str, dict):
        return params_str # already parsed
    if not isinstance(params_str, str):
        raise ValueError(f"params_str must be 'str', but received: {type(params_str)}")
    if params_str.startswith('['):
        params_str = params_str[1:]
    if params_str.endswith(']'):
        param_str = params_str[:-1]
    params = dict()
    for part in params_str.strip().split(','):
        k, v = part.split('='); params[k.strip()] = v.strip()
    return params

class Node:
    '''``Node`` class'''
    def __init__(self, label=None, edge_length=None):
        '''``Node`` constructor

        Args:
            ``label`` (``str``): Label of this ``Node``

            ``edge_length`` (``float``): Length of the edge incident to this ``Node``

        Returns:
            ``Node`` object
        '''
        self.children = []             # list of child Node objects
        self.parent = None             # parent Node object (None for root)
        self.label = label             # label
        self.edge_length = edge_length # length of incident edge

    def __lt__(self, other):
        '''Less Than operator. Just compares labels'''
        if not isinstance(other,Node):
            raise TypeError(f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'")
        elif self.label is None and other.label is not None:
            return True
        elif other.label is None:
            return False
        try:
            return float(self.label) < float(other.label)
        except:
            return str(self.label) < str(other.label)

    def __str__(self):
        '''Represent ``Node`` as a string (currently returns ``Node`` label as a string)

        Returns:
            ``str``: string representation of this ``Node``
        '''
        if self.label is None:
            return ''
        else:
            return str(self.label)

    def __copy__(self):
        '''Copy this ``Node``

        Returns:
            ``Node``: A copy of this ``Node``
        '''
        out = Node(label=copy(self.label), edge_length=copy(self.edge_length))
        out.children = copy(self.children)
        out.parent = self.parent
        return out

    def add_child(self, child):
        '''Add child to ``Node`` object

        Args:
            ``child`` (``Node``): The child ``Node`` to be added
        '''
        if not isinstance(child, Node):
            raise TypeError("child must be a Node")
        self.children.append(child); child.parent = self

    def child_nodes(self):
        '''Return a ``list`` containing this ``Node`` object's children

        Returns:
            ``list``: A ``list`` containing this ``Node`` object's children
        '''
        return copy(self.children)

    def contract(self):
        '''Contract this ``Node`` by directly connecting its children to its parent'''
        if self.is_root():
            return
        for c in self.children:
            if self.edge_length is not None and c.edge_length is not None:
                c.edge_length += self.edge_length
            self.parent.add_child(c)
        self.parent.remove_child(self)

    def get_edge_length(self):
        '''Return the length of the edge incident to this ``Node``

        Returns:
            ``float``: The length of the edge incident to this ``Node``
        '''
        return self.edge_length

    def get_label(self):
        '''Return the label of this ``Node``

        Returns:
            ``object``: The label of this ``Node``
        '''
        return self.label

    def get_parent(self):
        '''Return the parent of this ``Node``

        Returns:
            ``Node``: The parent of this ``Node``
        '''
        return self.parent

    def is_leaf(self):
        '''Returns ``True`` if this is a leaf

        Returns:
            ``bool``: ``True`` if this is a leaf, otherwise ``False``
        '''
        return len(self.children) == 0

    def is_root(self):
        '''Returns ``True`` if this is the ``root``

        Returns:
            ``bool``: ``True`` if this is the root, otherwise ``False``
        '''
        return self.parent is None

    def newick(self):
        '''Newick string conversion starting at this ``Node`` object

        Returns:
            ``str``: Newick string conversion starting at this ``Node`` object
        '''
        # set things up
        out = StringIO()
        root = self
        stack = [[root, 0]]

        # helper function to get node label as string
        def label_str(node):
            if node.label is None:
                return ''
            s = str(node.label)
            for c in UNSAFE_SYMBOLS:
                if c in s:
                    s = f"'{s}'"; break
            return s

        # helper function to append branch
        def append_branch(node):
            if hasattr(node, 'node_params'):
                out.write(params_str(node.node_params))
            if node.edge_length is not None or hasattr(node, 'edge_params'):
                out.write(':')
            if hasattr(node, 'edge_params'):
                out.write(params_str(node.edge_params))
            if isinstance(node.edge_length, float) and node.edge_length.is_integer():
                out.write(str(int(node.edge_length)))
            elif node.edge_length is not None:
                out.write(str(node.edge_length))
    
        # build Newick string
        while len(stack) != 0:
            node, i = stack[-1]
            children = node.children
            if i == 0:
                if len(children) == 0:
                    out.write(label_str(node))
                    if node is not root:
                        append_branch(node)
                    stack.pop(); continue
                out.write('(')
            if i < len(children):
                if i > 0:
                    out.write(',')
                child = children[i]
                stack[-1][1] += 1
                if len(child.children) == 0:
                    out.write(label_str(child))
                    append_branch(child)
                else:
                    stack.append([child, 0])
            else:
                out.write(')')
                if node.label is not None:
                    out.write(label_str(node))
                if node is not root:
                    append_branch(node)
                stack.pop()
        return out.getvalue()

    def num_children(self):
        '''Returns the number of children of this ``Node``

        Returns:
            ``int``: The number of children of this ``Node``
        '''
        return len(self.children)

    def num_nodes(self, leaves=True, internal=True):
        '''Compute the total number of selected nodes in the subtree rooted by this ``Node`` (including itself)

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

        Returns:
            ``int``: The total number of selected nodes in this ``Tree``
        '''
        if not isinstance(leaves, bool):
            raise TypeError("leaves must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        return sum((leaves and node.is_leaf()) or (internal and not node.is_leaf()) for node in self.traverse_preorder())

    def parse_params(self, node_params=True, edge_params=True):
        '''Attempt to parse the node and/or edge params of this ``Node`` as a ``dict``.

        Args:
            ``node_params`` (``bool``): ``True`` to attempt to parse this ``Node``'s ``node_params``

            ``edge_params`` (``bool``): ``True`` to attempt to parse this ``Node``'s ``edge_params``
        '''
        if node_params and hasattr(self, 'node_params') and isinstance(self.node_params, str):
            self.node_params = parse_params_str(self.node_params)
        if edge_params and hasattr(self, 'edge_params') and isinstance(self.edge_params, str):
            self.edge_params = parse_params_str(self.edge_params)

    def remove_child(self, child):
        '''Remove child from ``Node`` object

        Args:
            ``child`` (``Node``): The child to remove
        '''
        if not isinstance(child, Node):
            raise TypeError("child must be a Node")
        try:
            self.children.remove(child); child.parent = None
        except:
            raise RuntimeError("Attempting to remove non-existent child")

    def resolve_polytomies(self):
        '''Arbitrarily resolve polytomies below this ``Node`` with 0-lengthed edges.'''
        q = deque(); q.append(self)
        while len(q) != 0:
            node = q.popleft()
            while len(node.children) > 2:
                c1 = node.children.pop(); c2 = node.children.pop()
                nn = Node(edge_length=0); node.add_child(nn)
                nn.add_child(c1); nn.add_child(c2)
            q.extend(node.children)

    def set_edge_length(self, length):
        '''Set the length of the edge incident to this ``Node``

        Args:
            ``length``: The new length of the edge incident to this ``Node``
        '''
        try:
            self.edge_length = float(length)
        except:
            raise TypeError("length must be a float")

    def set_label(self, label):
        '''Set the label of this ``Node`` object

        Args:
            ``label``: The new label
        '''
        self.label = label

    def set_parent(self, parent):
        '''Set the parent of this ``Node`` object. Use this carefully, otherwise you may damage the structure of this ``Tree`` object.

        Args:
            ``Node``: The new parent of this ``Node``
        '''
        if not isinstance(parent, Node):
            raise TypeError("parent must be a Node")
        self.parent = parent

    def traverse_ancestors(self, include_self=True):
        '''Traverse over the ancestors of this ``Node``

        Args:
            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        if not isinstance(include_self, bool):
            raise TypeError("include_self must be a bool")
        if include_self:
            c = self
        else:
            c = self.parent
        while c is not None:
            yield c; c = c.parent

    def traverse_bfs(self, include_self=True):
        '''Perform a Breadth-First Search (BFS) starting at this ``Node`` object'. Yields (``Node``, distance) tuples
        
        Args:
            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        q = deque(); dist = {self: 0}; q.append((self,0))
        while len(q) != 0:
            curr = q.popleft()
            if include_self or curr[0] is not self:
                yield curr
            for c in curr[0].children:
                if c not in dist:
                    if c.edge_length is None:
                        el = 0
                    else:
                        el = c.edge_length
                    dist[c] = dist[curr[0]] + el; q.append((c,dist[c]))
            if curr[0].parent is not None and curr[0].parent not in dist:
                if curr[0].edge_length is None:
                    el = 0
                else:
                    el = curr[0].edge_length
                dist[curr[0].parent] = dist[curr[0]] + el; q.append((curr[0].parent,dist[curr[0].parent]))

    def traverse_inorder(self, leaves=True, internal=True, include_self=True):
        '''Perform an inorder traversal starting at this ``Node`` object

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        c = self; s = deque(); done = False
        while not done:
            if c is None:
                if len(s) == 0:
                    done = True
                else:
                    c = s.pop()
                    if ((leaves and c.is_leaf()) or (internal and not c.is_leaf())) and (include_self or c is not self):
                        yield c
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

    def traverse_internal(self, include_self=True):
        '''Traverse over the internal nodes below (and including) this ``Node`` object

        Args:
            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        yield from self.traverse_preorder(leaves=False, include_self=include_self)

    def traverse_leaves(self):
        '''Traverse over the leaves below this ``Node`` object

        Yields:
            ``Node``: The next node in the traversal
        '''
        yield from self.traverse_preorder(internal=False)

    def traverse_levelorder(self, leaves=True, internal=True, include_self=True):
        '''Perform a levelorder traversal starting at this ``Node`` object

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        q = deque(); q.append(self)
        while len(q) != 0:
            n = q.popleft()
            if ((leaves and n.is_leaf()) or (internal and not n.is_leaf())) and (include_self or n is not self):
                yield n
            q.extend(n.children)

    def traverse_postorder(self, leaves=True, internal=True, include_self=True):
        '''Perform a postorder traversal starting at this ``Node`` object

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        s1 = deque(); s2 = deque(); s1.append(self)
        while len(s1) != 0:
            n = s1.pop(); s2.append(n); s1.extend(n.children)
        while len(s2) != 0:
            n = s2.pop()
            if ((leaves and n.is_leaf()) or (internal and not n.is_leaf())) and (include_self or n is not self):
                yield n

    def traverse_preorder(self, leaves=True, internal=True, include_self=True):
        '''Perform a preorder traversal starting at this ``Node`` object

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``Node``: The next node in the traversal
        '''
        s = deque(); s.append(self)
        while len(s) != 0:
            n = s.pop()
            if ((leaves and n.is_leaf()) or (internal and not n.is_leaf())) and (include_self or n is not self):
                yield n
            s.extend(n.children)

    def traverse_rootdistorder(self, ascending=True, leaves=True, internal=True, include_self=True):
        '''Perform a traversal of the ``Node`` objects in the subtree rooted at this ``Node`` in either ascending (``ascending=True``) or descending (``ascending=False``) order of distance from this ``Node``

        Args:
            ``ascending`` (``bool``): ``True`` to perform traversal in ascending distance from the root, otherwise ``False`` for descending

            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``include_self`` (``bool``): ``True`` to include self in the traversal, otherwise ``False``

        Yields:
            ``tuple``: The next (root distance, ``Node``) pair in the traversal
        '''
        if not isinstance(ascending, bool):
            raise TypeError("ascending must be a bool")
        nodes = []; dist_from_root = {}
        for node in self.traverse_preorder():
            if node == self:
                d = 0
            else:
                d = dist_from_root[node.parent]
                if node.edge_length is not None:
                    d += node.edge_length
            dist_from_root[node] = d
            if ((leaves and node.is_leaf()) or (internal and not node.is_leaf())) and (include_self or node is not self):
                nodes.append((d,node))
        nodes.sort(reverse=(not ascending))
        yield from nodes
