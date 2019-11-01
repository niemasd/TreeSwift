from unittest import TestCase, main
import treeswift as ts


class TestNode(TestCase):

    def test_leaf_dijkstra(self):
        tree = ts.read_tree_newick("(A:3.2,(B:2.1,(C:1,D:1)));")
        nodes = [n for n in tree.traverse_preorder()]

        b = nodes[5]
        obs = b.leaf_dijkstra()
        exp = [(3.1,"C"), (3.1,"D"), (5.300000000000001,"A")]
        self.assertEqual(obs, exp)
        
        a = nodes[6]
        obs = a.leaf_dijkstra(2)
        exp = [(4.2,"C"), (4.2,"D")]
        self.assertEqual(obs, exp)
        
        tree = ts.read_tree_newick("((A:2.3,E:3)(B:2,(C:1.2,D:1)));")
        nodes = [n for n in tree.traverse_preorder()]

        a = nodes[7]
        obs = a.leaf_dijkstra(3)
        exp = [(3.3,"D"), (3.5,"C"), (4.3,"B")]
        self.assertEqual(obs, exp)

        with self.assertRaises(TypeError):
          obs = a.leaf_dijkstra(nodes[0])
        
        
      
 
if __name__ == '__main__':
    main()













