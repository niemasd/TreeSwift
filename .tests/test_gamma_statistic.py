"""Numerical regressions for Pybus--Harvey internode intervals.

Run with the local package on PYTHONPATH:
    python3 .tests/test_gamma_statistic.py
"""
from math import sqrt
import unittest
from treeswift import read_tree_newick


class GammaStatisticTests(unittest.TestCase):
    def assert_gamma(self, text, expected):
        tree = read_tree_newick(text)
        before = tree.newick()
        self.assertAlmostEqual(tree.gamma_statistic(), expected, places=12)
        self.assertEqual(tree.newick(), before)

    def test_known_interval_values(self):
        # For 3 tips, gamma = sqrt(3)*(2*g2-3*g3)/(2*g2+3*g3).
        cases = [
            ('((A:1,B:1)I:9,C:10)R;', 5*sqrt(3)/7),
            ('((A:1,B:1)I:1,C:2)R;', -sqrt(3)/5),
            ('((A:2,B:2)I:3,C:5)R;', 0),
            ('((A:9,B:9)I:1,C:10)R;', -25*sqrt(3)/29),
            ('((A:0,B:0)I:2,C:2)R;', sqrt(3)),
            # 4-tip intervals (1,1,1), (1,0,1), and (9,1,1).
            ('((A:1,B:1)I:2,(C:2,D:2)J:1)R;', -sqrt(24)/9),
            ('((A:1,B:1)I:1,(C:1,D:1)J:1)R;', -sqrt(24)/6),
            ('(((A:1,B:1)I:1,C:2)J:9,D:11)R;', 7*sqrt(24)/25),
        ]
        for text,expected in cases:
            with self.subTest(tree=text):
                self.assert_gamma(text, expected)

    def test_labels_do_not_change_branching_events(self):
        expected = 5*sqrt(3)/7
        for text in ['((A:1,B:1):9,C:10);', '((A:1,B:1)X:9,C:10);',
                     '((A:1,B:1):9,C:10)X;', '((:1,:1):9,:10);',
                     '((X:1,X:1)X:9,X:10)X;']:
            with self.subTest(tree=text):
                self.assert_gamma(text, expected)

    def test_scale_and_root_stem_invariance(self):
        for scale in [0.125, 1, 8]:
            for stem in [0, 17.25, 1e100]:
                text = '((A:%s,B:%s):%s,C:%s):%s;' % (
                    scale, scale, 9*scale, 10*scale, stem)
                with self.subTest(scale=scale, stem=stem):
                    self.assert_gamma(text, 5*sqrt(3)/7)

    def test_unary_nodes_are_not_speciation_events(self):
        # Prefix and edge subdivisions preserve the same two crown intervals.
        self.assert_gamma('(((A:0.5)U:0.5,B:1):9,C:10);', 5*sqrt(3)/7)
        self.assert_gamma('(((A:1,B:1):9,C:10):123)U:1e100;', 5*sqrt(3)/7)

    def test_polytomies_match_zero_length_binary_resolutions(self):
        for text in ['(A:3,B:3,C:3,D:3);',
                     '(A:3,(B:3,(C:3,D:3):0):0);']:
            self.assert_gamma(text, -sqrt(6))
        self.assert_gamma('(A:3,B:3,(C:1,D:1):2);', -sqrt(24)/5)

    def test_too_few_tips(self):
        for text in ['A;', '(A:1,B:1);', '(A:1,B:1)R:17;']:
            with self.subTest(tree=text):
                with self.assertRaisesRegex(RuntimeError, 'more than 2 leaves'):
                    read_tree_newick(text).gamma_statistic()

    def test_zero_total_remains_undefined(self):
        with self.assertRaises(ZeroDivisionError):
            read_tree_newick('((A:0,B:0)I:0,C:0)R;').gamma_statistic()


if __name__ == '__main__':
    unittest.main()
