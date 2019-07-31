#! /usr/bin/env python
from treeswift.Node import Node
from collections import deque
from copy import copy
from gzip import open as gopen
from math import ceil,log
from os.path import expanduser,isfile
from sys import version_info
from warnings import warn
INVALID_NEWICK = "Tree not valid Newick tree"
INVALID_NEXML = "Invalid NeXML file"
INVALID_NEXUS = "Invalid Nexus file"
EULER_GAMMA = 0.5772156649015328606065120900824024310421

# store bracket open/close for convenience in label parsing
BRACKET = {
    '[': ']', # square bracket
    '{': '}', # curly bracket
    "'": "'", # single-quote
    '"': '"', # double-quote
}

class Tree:
    '''``Tree`` class'''
    def __init__(self, is_rooted=True):
        '''``Tree`` constructor'''
        if not isinstance(is_rooted,bool):
            raise TypeError("is_rooted must be a bool")
        self.root = Node()  # root Node object
        self.is_rooted = is_rooted # boolean to see if the tree is rooted or not

    def __str__(self):
        '''Represent this ``Tree`` as a string

        Returns:
            ``str``: string representation of this ``Tree`` (Newick string)
        '''
        return self.newick()

    def __copy__(self):
        '''Copy this ``Tree``

        Returns:
            ``Tree``: A copy of this tree
        '''
        return self.extract_tree(None, False, False)


    def avg_branch_length(self, terminal=True, internal=True):
        '''Compute the average length of the selected branches of this ``Tree``. Edges with length ``None`` will be treated as 0-length

        Args:
            ``terminal`` (``bool``): ``True`` to include terminal branches, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal branches, otherwise ``False``

        Returns:
            The average length of the selected branches
        '''
        if not isinstance(terminal, bool):
            raise TypeError("terminal must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        if not internal and not terminal:
            raise RuntimeError("Must select either internal or terminal branches (or both)")
        tot = 0.; num = 0
        for node in self.traverse_preorder():
            if node.edge_length is not None and (internal and not node.is_leaf()) or (terminal and node.is_leaf()):
                tot += node.edge_length; num += 1
        return tot/num

    def branch_lengths(self, terminal=True, internal=True):
        '''Generator over the lengths of the selected branches of this ``Tree``. Edges with length ``None`` will be output as 0-length

        Args:
            ``terminal`` (``bool``): ``True`` to include terminal branches, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal branches, otherwise ``False``
        '''
        if not isinstance(terminal, bool):
            raise TypeError("terminal must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        for node in self.traverse_preorder():
            if (internal and not node.is_leaf()) or (terminal and node.is_leaf()):
                if node.edge_length is None:
                    yield 0
                else:
                    yield node.edge_length

    def closest_leaf_to_root(self):
        '''Return the leaf that is closest to the root and the corresponding distance. Edges with no length will be considered to have a length of 0

        Returns:
            ``tuple``: First value is the closest leaf to the root, and second value is the corresponding distance
        '''
        best = (None,float('inf')); d = dict()
        for node in self.traverse_preorder():
            if node.edge_length is None:
                d[node] = 0
            else:
                d[node] = node.edge_length
            if not node.is_root():
                d[node] += d[node.parent]
            if node.is_leaf() and d[node] < best[1]:
                best = (node,d[node])
        return best

    def coalescence_times(self, backward=True):
        '''Generator over the times of successive coalescence events

        Args:
            ``backward`` (``bool``): ``True`` to go backward in time (i.e., leaves to root), otherwise ``False``
        '''
        if not isinstance(backward, bool):
            raise TypeError("backward must be a bool")
        for dist in sorted((d for n,d in self.distances_from_root() if len(n.children) > 1), reverse=backward):
            yield dist

    def coalescence_waiting_times(self, backward=True):
        '''Generator over the waiting times of successive coalescence events

        Args:
            ``backward`` (``bool``): ``True`` to go backward in time (i.e., leaves to root), otherwise ``False``
        '''
        if not isinstance(backward, bool):
            raise TypeError("backward must be a bool")
        times = list(); lowest_leaf_dist = float('-inf')
        for n,d in self.distances_from_root():
            if len(n.children) > 1:
                times.append(d)
            elif len(n.children) == 0 and d > lowest_leaf_dist:
                lowest_leaf_dist = d
        times.append(lowest_leaf_dist)
        times.sort(reverse=backward)
        for i in range(len(times)-1):
            yield abs(times[i]-times[i+1])

    def collapse_short_branches(self, threshold):
        '''Collapse internal branches (not terminal branches) with length less than or equal to ``threshold``. A branch length of ``None`` is considered 0

        Args:
            ``threshold`` (``float``): The threshold to use when collapsing branches
        '''
        if not isinstance(threshold,float) and not isinstance(threshold,int):
            raise RuntimeError("threshold must be an integer or a float")
        elif threshold < 0:
            raise RuntimeError("threshold cannot be negative")
        q = deque(); q.append(self.root)
        while len(q) != 0:
            next = q.popleft()
            if next.edge_length is None or next.edge_length <= threshold:
                if next.is_root():
                    next.edge_length = None
                elif not next.is_leaf():
                    parent = next.parent; parent.remove_child(next)
                    for c in next.children:
                        parent.add_child(c)
            q.extend(next.children)

    def colless(self, normalize='leaves'):
        '''Compute the Colless balance index of this ``Tree``. If the tree has polytomies, they will be randomly resolved

        Args:
            ``normalize`` (``str``): How to normalize the Colless index (if at all)

            * ``None`` to not normalize

            * ``"leaves"`` to normalize by the number of leaves

            * ``"yule"`` to normalize to the Yule model

            * ``"pda"`` to normalize to the Proportional to Distinguishable Arrangements model

        Returns:
            ``float``: Colless index (either normalized or not)
        '''
        t_res = copy(self); t_res.resolve_polytomies(); leaves_below = dict(); n = 0; I = 0
        for node in t_res.traverse_postorder():
            if node.is_leaf():
                leaves_below[node] = 1; n += 1
            else:
                cl,cr = node.children; nl = leaves_below[cl]; nr = leaves_below[cr]
                leaves_below[node] = nl+nr; I += abs(nl-nr)
        if normalize is None or normalize is False:
            return I
        elif not isinstance(normalize,str):
            raise TypeError("normalize must be None or a string")
        normalize = normalize.lower()
        if normalize == 'leaves':
            return (2.*I)/((n-1)*(n-2))
        elif normalize == 'yule':
            return (I - n*log(n) - n*(EULER_GAMMA-1-log(2)))/n
        elif normalize == 'pda':
            return I/(n**1.5)
        else:
            raise RuntimeError("normalize must be None, 'leaves', 'yule', or 'pda'")

    def condense(self):
        '''If siblings have the same label, merge them. If they have edge lengths, the resulting ``Node`` will have the larger of the lengths'''
        self.resolve_polytomies(); labels_below = dict(); longest_leaf_dist = dict()
        for node in self.traverse_postorder():
            if node.is_leaf():
                labels_below[node] = [node.label]; longest_leaf_dist[node] = None
            else:
                labels_below[node] = set()
                for c in node.children:
                    labels_below[node].update(labels_below[c])
                    d = longest_leaf_dist[c]
                    if c.edge_length is not None:
                        if d is None:
                            d = 0
                        d += c.edge_length
                    if node not in longest_leaf_dist or longest_leaf_dist[node] is None or (d is not None and d > longest_leaf_dist[node]):
                        longest_leaf_dist[node] = d
        nodes = deque(); nodes.append(self.root)
        while len(nodes) != 0:
            node = nodes.pop()
            if node.is_leaf():
                continue
            if len(labels_below[node]) == 1:
                node.label = labels_below[node].pop(); node.children = list()
                if longest_leaf_dist[node] is not None:
                    if node.edge_length is None:
                        node.edge_length = 0
                    node.edge_length += longest_leaf_dist[node]
            else:
                nodes.extend(node.children)

    def contract_low_support(self, threshold):
        '''Contract internal nodes labeled by a number (e.g. branch support) below ``threshold``

        Args:
            ``threshold`` (``float``): The support threshold to use when contracting nodes'''
        if not isinstance(threshold, float) and not isinstance(threshold, int):
            raise TypeError("threshold must be float or int")
        to_contract = list()
        for node in self.traverse_preorder():
            try:
                if float(str(node)) < threshold:
                    to_contract.append(node)
            except:
                pass
        for node in to_contract:
            node.contract()

    def deroot(self, label='OLDROOT'):
        '''If the tree has a root edge, drop the edge to be a child of the root node

        Args:
            ``label`` (``str``): The desired label of the new child
        '''
        if self.root.edge_length is not None:
            self.root.add_child(Node(edge_length=self.root.edge_length,label=label))
            self.root.edge_length = None

    def diameter(self):
        '''Compute the diameter (maximum leaf pairwise distance) of this ``Tree``

        Returns:
            ``float``: The diameter of this Tree
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

    def distance_between(self, u, v):
        '''Return the distance between nodes ``u`` and ``v`` in this ``Tree``

        Args:
            ``u`` (``Node``): Node ``u``

            ``v`` (``Node``): Node ``v``

        Returns:
            ``float``: The distance between nodes ``u`` and ``v``
        '''
        if not isinstance(u, Node):
            raise TypeError("u must be a Node")
        if not isinstance(v, Node):
            raise TypeError("v must be a Node")
        if u == v:
            return 0.
        u_dists = {u:0.}; v_dists = {v:0.}
        c = u; p = u.parent # u traversal
        while p is not None:
            u_dists[p] = u_dists[c]
            if c.edge_length is not None:
                u_dists[p] += c.edge_length
            c = p; p = p.parent
        c = v; p = v.parent # v traversal
        while p is not None:
            v_dists[p] = v_dists[c]
            if c.edge_length is not None:
                v_dists[p] += c.edge_length
            if p in u_dists:
                return u_dists[p] + v_dists[p]
            c = p; p = p.parent
        raise RuntimeError("u and v are not in the same Tree")

    def distance_matrix(self, leaf_labels=False):
        '''Return a distance matrix (2D dictionary) of the leaves of this ``Tree``

        Args:
            ``leaf_labels`` (``bool``): ``True`` to have keys be labels of leaf ``Node`` objects, otherwise ``False`` to have keys be ``Node`` objects

        Returns:
            ``dict``: Distance matrix (2D dictionary) of the leaves of this ``Tree``, where keys are labels of leaves; ``M[u][v]`` = distance from ``u`` to ``v``
        '''
        M = dict(); leaf_dists = dict()
        for node in self.traverse_postorder():
            if node.is_leaf():
                leaf_dists[node] = [[node,0]]
            else:
                for c in node.children:
                    if c.edge_length is not None:
                        for i in range(len(leaf_dists[c])):
                            leaf_dists[c][i][1] += c.edge_length
                for c1 in range(0,len(node.children)-1):
                    leaves_c1 = leaf_dists[node.children[c1]]
                    for c2 in range(c1+1,len(node.children)):
                        leaves_c2 = leaf_dists[node.children[c2]]
                        for i in range(len(leaves_c1)):
                            for j in range(len(leaves_c2)):
                                u,ud = leaves_c1[i]; v,vd = leaves_c2[j]; d = ud+vd
                                if leaf_labels:
                                    u_key = u.label; v_key = v.label
                                else:
                                    u_key = u; v_key = v
                                if u_key not in M:
                                    M[u_key] = dict()
                                M[u_key][v_key] = d
                                if v_key not in M:
                                    M[v_key] = dict()
                                M[v_key][u_key] = d
                leaf_dists[node] = leaf_dists[node.children[0]]; del leaf_dists[node.children[0]]
                for i in range(1,len(node.children)):
                    leaf_dists[node] += leaf_dists[node.children[i]]; del leaf_dists[node.children[i]]
        return M

    def distances_from_parent(self, leaves=True, internal=True, unlabeled=False):
        '''Generator over the node-to-parent distances of this ``Tree``; (node,distance) tuples

        Args:
            ``terminal`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``unlabeled`` (``bool``): ``True`` to include unlabeled nodes, otherwise ``False``
        '''
        if not isinstance(leaves, bool):
            raise TypeError("leaves must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        if not isinstance(unlabeled, bool):
            raise TypeError("unlabeled must be a bool")
        if leaves or internal:
            for node in self.traverse_preorder():
                if ((leaves and node.is_leaf()) or (internal and not node.is_leaf())) and (unlabeled or node.label is not None):
                    if node.edge_length is None:
                        yield (node,0)
                    else:
                        yield (node,node.edge_length)

    def distances_from_root(self, leaves=True, internal=True, unlabeled=False):
        '''Generator over the root-to-node distances of this ``Tree``; (node,distance) tuples

        Args:
            ``terminal`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``

            ``unlabeled`` (``bool``): ``True`` to include unlabeled nodes, otherwise ``False``
        '''
        if not isinstance(leaves, bool):
            raise TypeError("leaves must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        if not isinstance(unlabeled, bool):
            raise TypeError("unlabeled must be a bool")
        if leaves or internal:
            d = dict()
            for node in self.traverse_preorder():
                if node.is_root():
                    d[node] = 0
                else:
                    d[node] = d[node.parent]
                if node.edge_length is not None:
                    d[node] += node.edge_length
                if ((leaves and node.is_leaf()) or (internal and not node.is_leaf())) and (unlabeled or node.label is not None):
                    yield (node,d[node])

    def edge_length_sum(self, terminal=True, internal=True):
        '''Compute the sum of all selected edge lengths in this ``Tree``

        Args:
            ``terminal`` (``bool``): ``True`` to include terminal branches, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal branches, otherwise ``False``

        Returns:
            ``float``: Sum of all selected edge lengths in this ``Tree``
        '''
        if not isinstance(terminal, bool):
            raise TypeError("leaves must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        return sum(node.edge_length for node in self.traverse_preorder() if node.edge_length is not None and ((terminal and node.is_leaf()) or (internal and not node.is_leaf())))

    def extract_subtree(self, node):
        '''Return a copy of the subtree rooted at ``node``

        Args:
            ``node`` (``Node``): The root of the desired subtree

        Returns:
            ``Tree``: A copy of the subtree rooted at ``node``
        '''
        if not isinstance(node, Node):
            raise TypeError("node must be a Node")
        r = self.root; self.root = node; o = copy(self); self.root = r; return o

    def extract_tree(self, labels, without, suppress_unifurcations=True):
        '''Helper function for ``extract_tree_*`` functions'''
        if not isinstance(suppress_unifurcations, bool):
            raise TypeError("suppress_unifurcations must be a bool")
        if labels is not None and not isinstance(labels, set):
            try:
                labels = set(labels)
            except:
                raise TypeError("labels must be iterable")
        label_to_leaf = dict(); keep = set()
        for node in self.traverse_leaves():
            label_to_leaf[str(node)] = node
            if labels is None or (without and str(node) not in labels) or (not without and str(node) in labels):
                keep.add(node)
        for node in list(keep):
            for a in node.traverse_ancestors(include_self=False):
                keep.add(a)
        out = Tree(); out.root.label = self.root.label; out.root.edge_length = self.root.edge_length
        q_old = deque(); q_old.append(self.root)
        q_new = deque(); q_new.append(out.root)
        while len(q_old) != 0:
            n_old = q_old.popleft(); n_new = q_new.popleft()
            for c_old in n_old.children:
                if c_old in keep:
                    c_new = Node(label=str(c_old), edge_length=c_old.edge_length); n_new.add_child(c_new)
                    q_old.append(c_old); q_new.append(c_new)
        if suppress_unifurcations:
            out.suppress_unifurcations()
        return out

    def extract_tree_without(self, labels, suppress_unifurcations=True):
        '''Extract a copy of this ``Tree`` without the leaves labeled by the strings in ``labels``

        Args:
            ``labels`` (``set``): Set of leaf labels to exclude

            ``suppress_unifurcations`` (``bool``): ``True`` to suppress unifurcations, otherwise ``False``

        Returns:
            ``Tree``: Copy of this ``Tree``, exluding the leaves labeled by the strings in ``labels``
        '''
        return self.extract_tree(labels, True, suppress_unifurcations)

    def extract_tree_with(self, labels, suppress_unifurcations=True):
        '''Extract a copy of this ``Tree`` with only the leaves labeled by the strings in ``labels``

        Args:
            ``leaves`` (``set``): Set of leaf labels to include.

            ``suppress_unifurcations`` (``bool``): ``True`` to suppress unifurcations, otherwise ``False``

        Returns:
            Tree: Copy of this Tree, including only the leaves labeled by the strings in ``labels``
        '''
        return self.extract_tree(labels, False, suppress_unifurcations)

    def furthest_from_root(self):
        '''Return the ``Node`` that is furthest from the root and the corresponding distance. Edges with no length will be considered to have a length of 0

        Returns:
            ``tuple``: First value is the furthest ``Node`` from the root, and second value is the corresponding distance
        '''
        best = (self.root,0); d = dict()
        for node in self.traverse_preorder():
            if node.edge_length is None:
                d[node] = 0
            else:
                d[node] = node.edge_length
            if not node.is_root():
                d[node] += d[node.parent]
            if d[node] > best[1]:
                best = (node,d[node])
        return best

    def gamma_statistic(self):
        '''Compute the Gamma statistic of Pybus and Harvey (2000)

        Returns:
            ``float``: The Gamma statistic of Pybus and Harvey (2000)
        '''
        t = copy(self); t.resolve_polytomies() # need fully bifurcating tree
        G = [g for g in t.coalescence_times(backward=False)]
        n = len(G)+1
        if n <= 2:
            raise RuntimeError("Gamma statistic can only be computed on trees with more than 2 leaves")
        T = sum((j+2)*g for j,g in enumerate(G))
        out = 0.
        for i in range(len(G)-1):
            for k in range(i+1):
                out += (k+2)*G[k]
        out /= (n-2)
        out -= (T/2)
        out /= T
        out /= (1./(12*(n-2)))**0.5
        return out

    def height(self):
        '''Compute the height (i.e., maximum distance from root) of this ``Tree``

        Returns:
            ``float``: The height (i.e., maximum distance from root) of this ``Tree``
        '''
        return max(d[1] for d in self.distances_from_root())

    def indent(self, space=4):
        '''Return an indented Newick string, just like ``nw_indent`` in Newick Utilities

        Args:
            ``space`` (``int``): The number of spaces a tab should equal

        Returns:
            ``str``: An indented Newick string
        '''
        if not isinstance(space,int):
            raise TypeError("space must be an int")
        if space < 0:
            raise ValueError("space must be a non-negative integer")
        space = ' '*space; o = []; l = 0
        for c in self.newick():
            if c == '(':
                o.append('(\n'); l += 1; o.append(space*l)
            elif c == ')':
                o.append('\n'); l -= 1; o.append(space*l); o.append(')')
            elif c == ',':
                o.append(',\n'); o.append(space*l)
            else:
                o.append(c)
        return ''.join(o)

    def label_to_node(self, selection='leaves'):
        '''Return a dictionary mapping labels (strings) to ``Node`` objects

        * If ``selection`` is ``"all"``, the dictionary will contain all nodes

        * If ``selection`` is ``"leaves"``, the dictionary will only contain leaves

        * If ``selection`` is ``"internal"``, the dictionary will only contain internal nodes

        * If ``selection`` is a ``set``, the dictionary will contain all nodes labeled by a label in ``selection``

        * If multiple nodes are labeled by a given label, only the last (preorder traversal) will be obtained

        Args:
            ``selection`` (``str`` or ``set``): The selection of nodes to get

            * ``"all"`` to select all nodes

            * ``"leaves"`` to select leaves

            * ``"internal"`` to select internal nodes

            * A ``set`` of labels to specify nodes to select

        Returns:
            ``dict``: Dictionary mapping labels to the corresponding nodes
        '''
        if not isinstance(selection,set) and not isinstance(selection,list) and (not isinstance(selection,str) or not (selection != 'all' or selection != 'leaves' or selection != 'internal')):
            raise RuntimeError('"selection" must be one of the strings "all", "leaves", or "internal", or it must be a set containing Node labels')
        if isinstance(selection, str):
            selection = selection[0]
        elif isinstance(selection,list):
            selection = set(selection)
        label_to_node = dict()
        for node in self.traverse_preorder():
            if selection == 'a' or (selection == 'i' and not node.is_leaf()) or (selection == 'l' and node.is_leaf()) or str(node) in selection:
                label_to_node[str(node)] = node
        if not isinstance(selection,str) and len(label_to_node) != len(selection):
            warn("Not all given labels exist in the tree")
        return label_to_node

    def labels(self, leaves=True, internal=True):
        '''Generator over the (non-``None``) ``Node`` labels of this ``Tree``

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``
        '''
        if not isinstance(leaves, bool):
            raise TypeError("leaves must be a bool")
        if not isinstance(internal, bool):
            raise TypeError("internal must be a bool")
        for node in self.traverse_preorder():
            if node.label is not None and ((leaves and node.is_leaf()) or (internal and not node.is_leaf())):
                yield node.label

    def ladderize(self, ascending=True):
        '''Ladderize this ``Tree`` by sorting each ``Node``'s children by total number of descendants

        Args:
            ``ascending`` (``bool``): ``True`` to sort in ascending order of ``mode``, otherwise ``False``
        '''
        self.order('num_descendants_then_edge_length_then_label', ascending=ascending)

    def lineages_through_time(self, present_day=None, show_plot=True, color='#000000', xmin=None, xmax=None, ymin=None, ymax=None, title=None, xlabel=None, ylabel=None):
        '''Compute the number of lineages through time. If seaborn is installed, a plot is shown as well

        Args:
            ``present_day`` (``float``): The time of the furthest node from the root. If ``None``, the top of the tree will be placed at time 0

            ``show_plot`` (``bool``): ``True`` to show the plot, otherwise ``False`` to only return the dictionary. To plot multiple LTTs on the same figure, set ``show_plot`` to False for all but the last plot

            ``color`` (``str``): The color of the resulting plot

            ``title`` (``str``): The title of the resulting plot

            ``xmin`` (``float``): The minimum value of the horizontal axis in the resulting plot

            ``xmax`` (``float``): The maximum value of the horizontal axis in the resulting plot

            ``xlabel`` (``str``): The label of the horizontal axis in the resulting plot

            ``ymin`` (``float``): The minimum value of the vertical axis in the resulting plot

            ``ymax`` (``float``): The maximum value of the vertical axis in the resulting plot

            ``ylabel`` (``str``): The label of the vertical axis in the resulting plot

        Returns:
            ``dict``: A dictionary in which each ``(t,n)`` pair denotes the number of lineages ``n`` that existed at time ``t``
        '''
        if present_day is not None and not isinstance(present_day,int) and not isinstance(present_day,float):
            raise TypeError("present_day must be a float")
        time = dict()
        if self.root.edge_length is None:
            tmproot = self.root
        else:
            tmproot = Node(); tmproot.add_child(self.root)
        for node in tmproot.traverse_preorder():
            if node.is_root():
                time[node] = 0.
            else:
                time[node] = time[node.parent]
                if node.edge_length is not None:
                    time[node] += node.edge_length
        nodes = sorted((time[node],node) for node in time)
        lineages = {nodes[0][0]:0}
        for i in range(len(nodes)):
            if nodes[i][0] not in lineages:
                lineages[nodes[i][0]] = lineages[nodes[i-1][0]]
            if nodes[i][1].edge_length is not None:
                if nodes[i][1].edge_length >= 0:
                    lineages[nodes[i][0]] -= 1
                else:
                    lineages[nodes[i][0]] += 1
            for c in nodes[i][1].children:
                if c.edge_length >= 0:
                    lineages[nodes[i][0]] += 1
                else:
                    lineages[nodes[i][0]] -= 1
        if present_day is not None:
            shift = present_day - max(lineages.keys())
        else:
            shift = max(0,-min(lineages.keys()))
        if shift != 0:
            tmp = dict()
            for t in lineages:
                tmp[t+shift] = lineages[t]
            lineages = tmp
        if tmproot != self.root:
            self.root.parent = None
        try:
            plot_ltt(lineages, show_plot=show_plot, color=color, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, title=title, xlabel=xlabel, ylabel=ylabel)
        except Exception as e:
            warn("Unable to produce visualization (but dictionary will still be returned)"); print(e)
        return lineages
    ltt = lineages_through_time # shorthand alias

    def mrca(self, labels):
        '''Return the Node that is the MRCA of the nodes labeled by a label in ``labels``. If multiple nodes are labeled by a given label, only the last (preorder traversal) will be obtained

        Args:
            ``labels`` (``set``): Set of leaf labels

        Returns:
            ``Node``: The MRCA of the ``Node`` objects labeled by a label in ``labels``
        '''
        if not isinstance(labels,set):
            try:
                labels = set(labels)
            except:
                raise TypeError("labels must be iterable")
        l2n = self.label_to_node(labels)
        count = dict()
        for node in l2n.values():
            for a in node.traverse_ancestors():
                if a not in count:
                    count[a] = 0
                count[a] += 1
                if count[a] == len(l2n):
                    return a
        raise RuntimeError("There somehow does not exist an MRCA for the given labels")

    def mrca_matrix(self):
        '''Return a dictionary storing all pairwise MRCAs. ``M[u][v]`` = MRCA of nodes ``u`` and ``v``. Excludes ``M[u][u]`` because MRCA of node and itself is itself

        Returns:
            ``dict``: ``M[u][v]`` = MRCA of nodes ``u`` and ``v``
        '''
        M = dict()
        leaves_below = dict()
        for node in self.traverse_postorder():
            leaves_below[node] = list()
            if node.is_leaf():
                leaves_below[node].append(node); M[node] = dict()
            else:
                for i in range(len(node.children)-1):
                    for l1 in leaves_below[node.children[i]]:
                        leaves_below[node].append(l1)
                        for j in range(i+1, len(node.children)):
                            for l2 in leaves_below[node.children[j]]:
                                M[l1][l2] = node; M[l2][l1] = node
                if len(node.children) != 1:
                    for l2 in leaves_below[node.children[-1]]:
                        leaves_below[node].append(l2)
        return M

    def newick(self):
        '''Output this ``Tree`` as a Newick string

        Returns:
            ``str``: Newick string of this ``Tree``
        '''
        if self.root.edge_length is None:
            suffix = ';'
        elif isinstance(self.root.edge_length,int):
            suffix = ':%d;' % self.root.edge_length
        elif isinstance(self.root.edge_length,float) and self.root.edge_length.is_integer():
            suffix = ':%d;' % int(self.root.edge_length)
        else:
            suffix = ':%s;' % str(self.root.edge_length)
        if self.is_rooted:
            return '[&R] %s%s' % (self.root.newick(),suffix)
        else:
            return '%s%s' % (self.root.newick(),suffix)

    def num_lineages_at(self, distance):
        '''Returns the number of lineages of this ``Tree`` that exist ``distance`` away from the root

        Args:
            ``distance`` (``float``): The distance away from the root

        Returns:
            ``int``: The number of lineages that exist ``distance`` away from the root
        '''
        if not isinstance(distance, float) and not isinstance(distance, int):
            raise TypeError("distance must be an int or a float")
        if distance < 0:
            raise RuntimeError("distance cannot be negative")
        d = dict(); q = deque(); q.append(self.root); count = 0
        while len(q) != 0:
            node = q.popleft()
            if node.is_root():
                d[node] = 0
            else:
                d[node] = d[node.parent]
            if node.edge_length is not None:
                d[node] += node.edge_length
            if d[node] < distance:
                q.extend(node.children)
            elif node.parent is None or d[node.parent] < distance:
                count += 1
        return count

    def num_nodes(self, leaves=True, internal=True):
        '''Compute the total number of selected nodes in this ``Tree``

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
        num = 0
        for node in self.traverse_preorder():
            if (leaves and node.is_leaf()) or (internal and not node.is_leaf()):
                num += 1
        return num

    def order(self, mode, ascending=True):
        '''Order the children of the nodes in this ``Tree`` based on ``mode``

        Args:
            ``mode`` (``str``): How to order the children of the nodes of this ``Tree``

            * ``"edge_length"`` = order by incident edge length

            * ``"edge_length_then_label"`` = order by incident edge length, then by node label

            * ``"edge_length_then_label_then_num_descendants"`` = order by incident edge length, then by node label, then by number of descendants

            * ``"edge_length_then_num_descendants"`` = order by incident edge length, then by number of descendants

            * ``"edge_length_then_num_descendants_then_label"`` = order by incident edge length, then by number of descendants, then by node label

            * ``"label"`` = order by node label

            * ``"label_then_edge_length"`` = order by node label, then by incident edge length

            * ``"label_then_edge_length_then_num_descendants"`` = order by node label, then by incident edge length, then by number of descendants

            * ``"label_then_num_descendants"`` = order by node label, then by number of descendants

            * ``"label_then_num_descendants_then_edge_length"`` = order by node label, then by number of descendants, then by incident edge length

            * ``"num_descendants"`` = order by number of descendants

            * ``"num_descendants_then_label"`` = order by number of descendants, then by node label

            * ``"num_descendants_then_label_then_edge_length"`` = order by number of descendants, then by node label, then by incident edge length

            * ``"num_descendants_then_edge_length"`` = order by number of descendants, then by incident edge length

            * ``"num_descendants_then_edge_length_then_label"`` = order by number of descendants, then by incident edge length, then by node label

            ``ascending`` (``bool``): ``True`` to sort in ascending order of ``mode``, otherwise ``False``
        '''
        if not isinstance(mode, str):
            raise TypeError("mode must be a str")
        if not isinstance(ascending, bool):
            raise TypeError("ascending must be a bool")
        if 'num_descendants' in mode:
            num_descendants = dict()
            for node in self.traverse_postorder():
                if node.is_leaf():
                    num_descendants[node] = 0
                else:
                    num_descendants[node] = sum(num_descendants[c] for c in node.children) + len(node.children)
        if mode == 'edge_length':
            k = lambda node: (node.edge_length is not None, node.edge_length)
        elif mode == 'edge_length_then_label':
            k = lambda node: (node.edge_length is not None, node.edge_length, node.label is not None, node.label)
        elif mode == 'edge_length_then_label_then_num_descendants':
            k = lambda node: (node.edge_length is not None, node.edge_length, node.label is not None, node.label, num_descendants[node])
        elif mode == 'edge_length_then_num_descendants':
            k = lambda node: (node.edge_length is not None, node.edge_length, num_descendants[node])
        elif mode == 'edge_length_then_num_descendants_then_label':
            k = lambda node: (node.edge_length is not None, node.edge_length, num_descendants[node], node.label is not None, node.label)
        elif mode == 'label':
            k = lambda node: (node.label is not None, node.label)
        elif mode == 'label_then_edge_length':
            k = lambda node: (node.label is not None, node.label, node.edge_length is not None, node.edge_length)
        elif mode == 'label_then_edge_length_then_num_descendants':
            k = lambda node: (node.label is not None, node.label, node.edge_length is not None, node.edge_length, num_descendants[node])
        elif mode == 'label_then_num_descendants':
            k = lambda node: (node.label is not None, node.label, num_descendants[node])
        elif mode == 'label_then_num_descendants_then_edge_length':
            k = lambda node: (node.label is not None, node.label, num_descendants[node], node.edge_length is not None, node.edge_length)
        elif mode == 'num_descendants':
            k = lambda node: num_descendants[node]
        elif mode == 'num_descendants_then_label':
            k = lambda node: (num_descendants[node], node.label is not None, node.label)
        elif mode == 'num_descendants_then_label_then_edge_length':
            k = lambda node: (num_descendants[node], node.label is not None, node.label, node.edge_length is not None, node.edge_length)
        elif mode == 'num_descendants_then_edge_length':
            k = lambda node: (num_descendants[node], node.edge_length is not None, node.edge_length)
        elif mode == 'num_descendants_then_edge_length_then_label':
            k = lambda node: (num_descendants[node], node.edge_length is not None, node.edge_length, node.label is not None, node.label)
        else:
            raise ValueError("Invalid choice for mode")
        for node in self.traverse_preorder():
            node.children.sort(key=k, reverse=not ascending)

    def rename_nodes(self, renaming_map):
        '''Rename nodes in this ``Tree``

        Args:
            ``renaming_map`` (``dict``): A dictionary mapping old labels (keys) to new labels (values)
        '''
        if not isinstance(renaming_map, dict):
            raise TypeError("renaming_map must be a dict")
        for node in self.traverse_preorder():
            if node.label in renaming_map:
                node.label = renaming_map[node.label]

    def reroot(self, node, length=None, branch_support=False):
        '''Reroot this ``Tree`` at ``length`` up the incident edge of ``node``. If 0 or ``None``, reroot at the node (not on the incident edge)

        Args:
            ``node`` (``Node``): The ``Node`` on whose incident edge this ``Tree`` will be rerooted

            ``length`` (``float``): The distance up the specified edge at which to reroot this ``Tree``. If 0 or ``None``, reroot at the node (not on the incident edge)

            ``branch_support`` (``bool``): ``True`` if internal node labels represent branch support values, otherwise ``False``
        '''
        warn("TreeSwift's rerooting functionality is poorly tested and likely has bugs. It will be fixed in a future release")
        if not isinstance(node, Node):
            raise TypeError("node must be a Node")
        if length is not None and not isinstance(length, float) and not isinstance(length, int):
            raise TypeError("length must be a float")
        if not isinstance(branch_support, bool):
            raise TypeError("branch_support must be a bool")
        if length is not None and length < 0:
            raise ValueError("Specified length at which to reroot must be positive")
        if node.edge_length is None:
            if length is not None and length != 0:
                raise ValueError("Specified node has no edge length, so specified length must be None or 0")
        elif length is not None and length > node.edge_length:
            raise ValueError("Specified length must be shorter than the edge at which to reroot")
        if length is not None and length > 0:
            newnode = Node(edge_length=node.edge_length-length); node.edge_length -= length
            if not node.is_root():
                p = node.parent; p.children.remove(node); p.add_child(newnode)
            newnode.add_child(node); node = newnode
        if node.is_root():
            return
        elif self.root.edge_length is not None:
            newnode = Node(label='ROOT'); newnode.add_child(self.root); self.root = newnode
        ancestors = [a for a in node.traverse_ancestors(include_self=True) if not a.is_root()]
        for i in range(len(ancestors)-1, -1, -1):
            curr = ancestors[i]; curr.parent.edge_length = curr.edge_length; curr.edge_length = None
            if branch_support:
                curr.parent.label = curr.label; curr.label = None
            curr.parent.children.remove(curr); curr.add_child(curr.parent); curr.parent = None
        self.root = node; self.is_rooted = True

    def resolve_polytomies(self):
        '''Arbitrarily resolve polytomies with 0-lengthed edges.'''
        self.root.resolve_polytomies()

    def sackin(self, normalize='leaves'):
        '''Compute the Sackin balance index of this ``Tree``

        Args:
            ``normalize`` (``str``): How to normalize the Sackin index (if at all)

            * ``None`` to not normalize

            * ``"leaves"`` to normalize by the number of leaves

            * ``"yule"`` to normalize to the Yule model

            * ``"pda"`` to normalize to the Proportional to Distinguishable Arrangements model

        Returns:
            ``float``: Sackin index (either normalized or not)
        '''
        num_nodes_from_root = dict(); sackin = 0; num_leaves = 0
        for node in self.traverse_preorder():
            num_nodes_from_root[node] = 1
            if not node.is_root():
                num_nodes_from_root[node] += num_nodes_from_root[node.parent]
            if node.is_leaf():
                num_nodes_from_root[node] -= 1; sackin += num_nodes_from_root[node]; num_leaves += 1
        if normalize is None or normalize is False:
            return sackin
        elif not isinstance(normalize,str):
            raise TypeError("normalize must be None or a string")
        normalize = normalize.lower()
        if normalize == 'leaves':
            return float(sackin)/num_leaves
        elif normalize == 'yule':
            x = sum(1./i for i in range(2, num_leaves+1))
            return (sackin - (2*num_leaves*x)) / num_leaves
        elif normalize == 'pda':
            return sackin/(num_leaves**1.5)
        else:
            raise RuntimeError("normalize must be None, 'leaves', 'yule', or 'pda'")

    def scale_edges(self, multiplier):
        '''Multiply all edges in this ``Tree`` by ``multiplier``'''
        if not isinstance(multiplier,int) and not isinstance(multiplier,float):
            raise TypeError("multiplier must be an int or float")
        for node in self.traverse_preorder():
            if node.edge_length is not None:
                node.edge_length *= multiplier

    def suppress_unifurcations(self):
        '''Remove all nodes with only one child and directly attach child to parent'''
        q = deque(); q.append(self.root)
        while len(q) != 0:
            node = q.popleft()
            if len(node.children) != 1:
                q.extend(node.children); continue
            child = node.children.pop()
            if node.is_root():
                self.root = child; child.parent = None
            else:
                parent = node.parent; parent.remove_child(node); parent.add_child(child)
            if node.edge_length is not None:
                if child.edge_length is None:
                    child.edge_length = 0
                child.edge_length += node.edge_length
            if child.label is None and node.label is not None:
                child.label = node.label
            q.append(child)

    def traverse_inorder(self, leaves=True, internal=True):
        '''Perform an inorder traversal of the ``Node`` objects in this ``Tree``

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``
        '''
        for node in self.root.traverse_inorder(leaves=leaves, internal=internal):
            yield node

    def traverse_internal(self):
        '''Traverse over the internal nodes of this ``Tree``'''
        for node in self.root.traverse_internal():
            yield node

    def traverse_leaves(self):
        '''Traverse over the leaves of this ``Tree``'''
        for node in self.root.traverse_leaves():
            yield node

    def traverse_levelorder(self, leaves=True, internal=True):
        '''Perform a levelorder traversal of the ``Node`` objects in this ``Tree``'''
        for node in self.root.traverse_levelorder(leaves=leaves, internal=internal):
            yield node

    def traverse_postorder(self, leaves=True, internal=True):
        '''Perform a postorder traversal of the ``Node`` objects in this ``Tree``

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``
        '''
        for node in self.root.traverse_postorder(leaves=leaves, internal=internal):
            yield node

    def traverse_preorder(self, leaves=True, internal=True):
        '''Perform a preorder traversal of the ``Node`` objects in this ``Tree``

        Args:
            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``
        '''
        for node in self.root.traverse_preorder(leaves=leaves, internal=internal):
            yield node

    def traverse_rootdistorder(self, ascending=True, leaves=True, internal=True):
        '''Perform a traversal of the ``Node`` objects in this ``Tree`` in either ascending (``ascending=True``) or descending (``ascending=False``) order of distance from the root

        Args:
            ``ascending`` (``bool``): ``True`` to perform traversal in ascending distance from the root, otherwise ``False`` for descending

            ``leaves`` (``bool``): ``True`` to include leaves, otherwise ``False``

            ``internal`` (``bool``): ``True`` to include internal nodes, otherwise ``False``
        '''
        for node in self.root.traverse_rootdistorder(ascending=ascending, leaves=leaves, internal=internal):
            yield node

    def treeness(self):
        '''Compute the `treeness` (sum of internal branch lengths / sum of all branch lengths) of this ``Tree``. Branch lengths of ``None`` are considered 0 length

        Returns:
            ``float``: `Treeness` of this ``Tree`` (sum of internal branch lengths / sum of all branch lengths)
        '''
        internal = 0.; all = 0.
        for node in self.traverse_preorder():
            if node.edge_length is not None:
                all += node.edge_length
                if not node.is_leaf():
                    internal += node.edge_length
        return internal/all

    def write_tree_newick(self, filename, hide_rooted_prefix=False):
        '''Write this ``Tree`` to a Newick file

        Args:
            ``filename`` (``str``): Path to desired output file (plain-text or gzipped)
        '''
        if not isinstance(filename, str):
            raise TypeError("filename must be a str")
        treestr = self.newick()
        if hide_rooted_prefix:
            if treestr.startswith('[&R]'):
                treestr = treestr[4:].strip()
            else:
                warn("Specified hide_rooted_prefix, but tree was not rooted")
        if filename.lower().endswith('.gz'): # gzipped file
            f = gopen(expanduser(filename),'wb',9); f.write(treestr.encode()); f.close()
        else: # plain-text file
            f = open(expanduser(filename),'w'); f.write(treestr); f.close()

def plot_ltt(lineages, show_plot=True, color='#000000', xmin=None, xmax=None, ymin=None, ymax=None, title=None, xlabel=None, ylabel=None):
    '''Plot the Lineages Through Time (LTT) curve of a given tree

    Args:
        ``lineages`` (``dict``): The ``lineages`` dictionary returned by a ``Tree`` object's ``lineages_through_time()`` function call

        ``show_plot`` (``bool``): ``True`` to show the plot, otherwise ``False`` to only return the dictionary. To plot multiple LTTs on the same figure, set ``show_plot`` to False for all but the last plot.

        ``color`` (``str``): The color of the resulting plot

        ``title`` (``str``): The title of the resulting plot

        ``xmin`` (``float``): The minimum value of the horizontal axis in the resulting plot

        ``xmax`` (``float``): The maximum value of the horizontal axis in the resulting plot

        ``xlabel`` (``str``): The label of the horizontal axis in the resulting plot

        ``ymin`` (``float``): The minimum value of the vertical axis in the resulting plot

        ``ymax`` (``float``): The maximum value of the vertical axis in the resulting plot

        ``ylabel`` (``str``): The label of the vertical axis in the resulting plot
    '''
    import matplotlib.pyplot as plt; from matplotlib.ticker import MaxNLocator
    if 'TREESWIFT_FIGURE' not in globals():
        global TREESWIFT_FIGURE; TREESWIFT_FIGURE = None
    if TREESWIFT_FIGURE is None:
        TREESWIFT_FIGURE = plt.figure()
        TREESWIFT_FIGURE.gca().yaxis.set_major_locator(MaxNLocator(integer=True)) # integer y ticks
        TREESWIFT_FIGURE.XMIN = float('inf'); TREESWIFT_FIGURE.XMAX = float('-inf')
        TREESWIFT_FIGURE.YMIN = float('inf'); TREESWIFT_FIGURE.YMAX = float('-inf')
    times = sorted(lineages.keys())
    if times[0] < TREESWIFT_FIGURE.XMIN:
        TREESWIFT_FIGURE.XMIN = times[0]
    if times[-1] > TREESWIFT_FIGURE.XMAX:
        TREESWIFT_FIGURE.XMAX = times[-1]
    for i in range(len(times)-1):
        if i == 0:
            prev = 0
        else:
            prev = lineages[times[i-1]]
        if lineages[times[i]] > TREESWIFT_FIGURE.YMAX:
            TREESWIFT_FIGURE.YMAX = lineages[times[i]]
        if lineages[times[i]] < TREESWIFT_FIGURE.YMIN:
            TREESWIFT_FIGURE.YMIN = lineages[times[i]]
        TREESWIFT_FIGURE.gca().plot([times[i],times[i]], [prev,lineages[times[i]]], color=color)
        TREESWIFT_FIGURE.gca().plot([times[i],times[i+1]], [lineages[times[i]],lineages[times[i]]], color=color)
    if len(times) > 1:
        TREESWIFT_FIGURE.gca().plot([times[-1],times[-1]], [lineages[times[-2]],lineages[times[-1]]], color=color)
        if lineages[times[-1]] < TREESWIFT_FIGURE.YMIN:
            TREESWIFT_FIGURE.YMIN = lineages[times[-1]]
    if show_plot:
        if xmin is None:
            xmin = TREESWIFT_FIGURE.XMIN
        elif not isinstance(xmin,int) and not isinstance(xmin,float):
            warn("xmin is invalid, so using the default"); xmin = TREESWIFT_FIGURE.XMIN
        if xmax is None:
            xmax = TREESWIFT_FIGURE.XMAX
        elif not isinstance(xmax,int) and not isinstance(xmax,float):
            warn("xmax is invalid, so using the default"); xmax = TREESWIFT_FIGURE.XMAX
        plt.xlim(left=xmin, right=xmax)
        if ymin is None:
            ymin = TREESWIFT_FIGURE.YMIN
        elif not isinstance(ymin,int) and not isinstance(ymin,float):
            warn("ymin is invalid, so using the default"); ymin = TREESWIFT_FIGURE.YMIN
        if ymax is None:
            ymax = ceil(TREESWIFT_FIGURE.YMAX*1.1)
        elif not isinstance(ymax,int) and not isinstance(ymax,float):
            warn("ymax is invalid, so using the default"); ymax = ceil(TREESWIFT_FIGURE.YMAX*1.1)
        plt.ylim(bottom=ymin, top=ymax)
        if title is not None and not isinstance(title,str):
            warn("title is invalid, so using the default"); title = None
        if title is None:
            plt.title("Lineages Through Time")
        else:
            plt.title(title)
        if xlabel is not None and not isinstance(xlabel,str):
            warn("xlabel is invalid, so using the default"); xlabel = None
        if xlabel is None:
            plt.xlabel("Time")
        else:
            plt.xlabel(xlabel)
        if ylabel is not None and not isinstance(ylabel,str):
            warn("ylabel is invalid, so using the default"); ylabel = None
        if ylabel is None:
            plt.ylabel("Number of Lineages")
        else:
            plt.ylabel(ylabel)
        plt.show(); TREESWIFT_FIGURE = None

def read_tree_dendropy(tree):
    '''Create a TreeSwift tree from a DendroPy tree

    Args:
        ``tree`` (``dendropy.datamodel.treemodel``): A Dendropy ``Tree`` object

    Returns:
        ``Tree``: A TreeSwift tree created from ``tree``
    '''
    out = Tree(); d2t = dict()
    if not hasattr(tree, 'preorder_node_iter') or not hasattr(tree, 'seed_node') or not hasattr(tree, 'is_rooted'):
        raise TypeError("tree must be a DendroPy Tree object")
    if tree.is_rooted != True:
        out.is_rooted = False
    for node in tree.preorder_node_iter():
        if node == tree.seed_node:
            curr = out.root
        else:
            curr = Node(); d2t[node.parent_node].add_child(curr)
        d2t[node] = curr; curr.edge_length = node.edge_length
        if hasattr(node, 'taxon') and node.taxon is not None:
            curr.label = node.taxon.label
        else:
            curr.label = node.label
    return out

def read_tree_newick(newick):
    '''Read a tree from a Newick string or file

    Args:
        ``newick`` (``str``): Either a Newick string or the path to a Newick file (plain-text or gzipped)

    Returns:
        ``Tree``: The tree represented by ``newick``. If the Newick file has multiple trees (one per line), a ``list`` of ``Tree`` objects will be returned
    '''
    if not isinstance(newick, str):
        try:
            newick = str(newick)
        except:
            raise TypeError("newick must be a str")
    if newick.lower().endswith('.gz'): # gzipped file
        f = gopen(expanduser(newick)); ts = f.read().decode().strip(); f.close()
    elif isfile(expanduser(newick)): # plain-text file
        f = open(expanduser(newick)); ts = f.read().strip(); f.close()
    else:
        ts = newick.strip()
    lines = ts.splitlines()
    if len(lines) != 1:
        return [read_tree_newick(l) for l in lines]
    try:
        t = Tree(); t.is_rooted = ts.startswith('[&R]')
        if ts[0] == '[':
            ts = ']'.join(ts.split(']')[1:]).strip(); ts = ts.replace(', ',',')
        n = t.root; i = 0
        while i < len(ts):
            # end of Newick string
            if ts[i] == ';':
                if i != len(ts)-1 or n != t.root:
                    raise RuntimeError(INVALID_NEWICK)

            # go to new child
            elif ts[i] == '(':
                c = Node(); n.add_child(c); n = c

            # go to parent
            elif ts[i] == ')':
                n = n.parent

            # go to new sibling
            elif ts[i] == ',':
                n = n.parent; c = Node(); n.add_child(c); n = c

            # edge length
            elif ts[i] == ':':
                i += 1; ls = ''
                while ts[i] != ',' and ts[i] != ')' and ts[i] != ';':
                    ls += ts[i]; i += 1
                if ls[0] == '[':
                    n.edge_params = ']'.join(ls.split(']')[:-1]); ls = ls.split(']')[-1]
                n.edge_length = float(ls); i -= 1

            # node label
            else:
                label = ''; bracket = None
                while bracket is not None or ts[i] in BRACKET or (ts[i] != ':' and ts[i] != ',' and ts[i] != ';' and ts[i] != ')'):
                    if ts[i] in BRACKET and bracket is None:
                        bracket = ts[i]
                    elif bracket is not None and ts[i] == BRACKET[bracket]:
                        bracket = None
                    label += ts[i]; i += 1
                i -= 1; n.label = label
            i += 1
    except Exception as e:
        raise RuntimeError("Failed to parse string as Newick: %s"%ts)
    return t

def read_tree_nexml(nexml):
    '''Read a tree from a NeXML string or file

    Args:
        ``nexml`` (``str``): Either a NeXML string or the path to a NeXML file (plain-text or gzipped)

    Returns:
        ``dict`` of ``Tree``: A dictionary of the trees represented by ``nexml``, where keys are tree names (``str``) and values are ``Tree`` objects
    '''
    if not isinstance(nexml, str):
        raise TypeError("nexml must be a str")
    if nexml.lower().endswith('.gz'): # gzipped file
        f = gopen(expanduser(nexml))
    elif isfile(expanduser(nexml)): # plain-text file
        f = open(expanduser(nexml))
    else:
        f = nexml.splitlines()
    trees = dict(); id_to_node = dict(); tree_id = None
    for line in f:
        if isinstance(line,bytes):
            l = line.decode().strip()
        else:
            l = line.strip()
        l_lower = l.lower()
        # start of tree
        if l_lower.startswith('<tree '):
            if tree_id is not None:
                raise ValueError(INVALID_NEXML)
            parts = l.split()
            for part in parts:
                if '=' in part:
                    k,v = part.split('='); k = k.strip()
                    if k.lower() == 'id':
                        tree_id = v.split('"')[1]; break
            if tree_id is None:
                raise ValueError(INVALID_NEXML)
            trees[tree_id] = Tree(); trees[tree_id].root = None
        # end of tree
        elif l_lower.replace(' ','').startswith('</tree>'):
            if tree_id is None:
                raise ValueError(INVALID_NEXML)
            id_to_node = dict(); tree_id = None
        # node
        elif l_lower.startswith('<node '):
            if tree_id is None:
                raise ValueError(INVALID_NEXML)
            node_id = None; node_label = None; is_root = False
            k = ''; v = ''; in_key = True; in_quote = False
            for i in range(6, len(l)):
                if l[i] == '"' or l[i] == "'":
                    in_quote = not in_quote
                if not in_quote and in_key and l[i] == '=':
                    in_key = False
                elif not in_quote and not in_key and (l[i] == '"' or l[i] == "'"):
                    k = k.strip()
                    if k.lower() == 'id':
                        node_id = v
                    elif k.lower() == 'label':
                        node_label = v
                    elif k.lower() == 'root' and v.strip().lower() == 'true':
                        is_root = True
                    in_key = True; k = ''; v = ''
                elif in_key and not (l[i] == '"' or l[i] == "'"):
                    k += l[i]
                elif not in_key and not (l[i] == '"' or l[i] == "'"):
                    v += l[i]
            if node_id is None or node_id in id_to_node:
                raise ValueError(INVALID_NEXML)
            id_to_node[node_id] = Node(label=node_label)
            if is_root:
                if trees[tree_id].root is not None:
                    raise ValueError(INVALID_NEXML)
                trees[tree_id].root = id_to_node[node_id]
        # edge
        elif l_lower.startswith('<edge '):
            if tree_id is None:
                raise ValueError(INVALID_NEXML)
            source = None; target = None; length = None
            parts = l.split()
            for part in parts:
                if '=' in part:
                    k,v = part.split('='); k = k.strip(); k_lower = k.lower()
                    if k_lower == 'source':
                        source = v.split('"')[1]
                    elif k_lower == 'target':
                        target = v.split('"')[1]
                    elif k_lower == 'length':
                        length = float(v.split('"')[1])
            if source is None or target is None or length is None:
                raise ValueError(INVALID_NEXML)
            if source not in id_to_node:
                raise ValueError(INVALID_NEXML)
            if target not in id_to_node:
                raise ValueError(INVALID_NEXML)
            id_to_node[source].add_child(id_to_node[target])
            id_to_node[target].edge_length = length
        elif l_lower.startswith('<rootedge '):
            if tree_id is None:
                raise ValueError(INVALID_NEXML)
            root_node = None; length = None
            parts = l.split()
            for part in parts:
                if '=' in part:
                    k,v = part.split('='); k = k.strip(); k_lower = k.lower()
                    if k_lower == 'target':
                        root_node = id_to_node[v.split('"')[1]]
                    elif k_lower == 'length':
                        length = float(v.split('"')[1])
            if trees[tree_id].root is None:
                raise ValueError(INVALID_NEXML)
            if root_node is not None and trees[tree_id].root != root_node:
                raise ValueError(INVALID_NEXML)
            trees[tree_id].root.edge_length = length
    if hasattr(f,'close'):
        f.close()
    return trees

def read_tree_nexus(nexus):
    '''Read a tree from a Nexus string or file

    Args:
        ``nexus`` (``str``): Either a Nexus string or the path to a Nexus file (plain-text or gzipped)

    Returns:
        ``dict`` of ``Tree``: A dictionary of the trees represented by ``nexus``, where keys are tree names (``str``) and values are ``Tree`` objects
    '''
    if not isinstance(nexus, str):
        raise TypeError("nexus must be a str")
    if nexus.lower().endswith('.gz'): # gzipped file
        f = gopen(expanduser(nexus))
    elif isfile(expanduser(nexus)): # plain-text file
        f = open(expanduser(nexus))
    else:
        f = nexus.splitlines()
    trees = dict()
    for line in f:
        if isinstance(line,bytes):
            l = line.decode().strip()
        else:
            l = line.strip()
        if l.lower().startswith('tree '):
            i = l.index('='); left = l[:i].strip(); right = l[i+1:].strip()
            name = ' '.join(left.split(' ')[1:])
            trees[name] = read_tree_newick(right)
    if hasattr(f,'close'):
        f.close()
    if len(trees) == 0:
        raise ValueError(INVALID_NEXUS)
    return trees

def read_tree(input, schema):
    '''Read a tree from a string or file

    Args:
        ``input`` (``str``): Either a tree string, a path to a tree file (plain-text or gzipped), or a DendroPy Tree object

        ``schema`` (``str``): The schema of ``input`` (DendroPy, Newick, NeXML, or Nexus)

    Returns:
        * If the input is Newick, either a ``Tree`` object if ``input`` contains a single tree, or a ``list`` of ``Tree`` objects if ``input`` contains multiple trees (one per line)

        * If the input is NeXML or Nexus, a ``dict`` of trees represented by ``input``, where keys are tree names (``str``) and values are ``Tree`` objects
    '''
    schema_to_function = {
        'dendropy': read_tree_dendropy,
        'newick': read_tree_newick,
        'nexml': read_tree_nexml,
        'nexus': read_tree_nexus
    }
    if schema.lower() not in schema_to_function:
        raise ValueError("Invalid schema: %s (valid options: %s)" % (schema, ', '.join(sorted(schema_to_function.keys()))))
    return schema_to_function[schema.lower()](input)
