# Copyright (c) 2020, CNRS
# Authors: Pierre Fernbach <pfernbac@laas.fr>

import pickle
import unittest

from numpy import array

from ndcurves import curve_constraints


class CurveConstraintsTest(unittest.TestCase):
    def test_operator_equal(self):
        c = curve_constraints(3)
        c.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_vel = array([[0.0, -1.0, 1.0]]).transpose()
        c.init_acc = array([[0.0, 1.0, -1.0]]).transpose()
        c.end_acc = array([[0.0, 100.0, 1.0]]).transpose()
        c.init_jerk = array([[2.0, 4.0, 1.0]]).transpose()
        c.end_jerk = array([[-1.0, 2.0, 7.0]]).transpose()

        c2 = curve_constraints(3)
        c2.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c2.end_vel = array([[0.0, -1.0, 1.0]]).transpose()
        c2.init_acc = array([[0.0, 1.0, -1.0]]).transpose()
        c2.end_acc = array([[0.0, 100.0, 1.0]]).transpose()
        c2.init_jerk = array([[2.0, 4.0, 1.0]]).transpose()
        c2.end_jerk = array([[-1.0, 2.0, 7.0]]).transpose()
        self.assertTrue(c == c2)
        c2.init_vel = array([[1.0, 1.0, 1.0]]).transpose()
        self.assertTrue(c != c2)

    def test_serialization(self):
        c = curve_constraints(3)
        c.init_vel = array([[0.0, 1.0, 1.0]]).transpose()
        c.end_vel = array([[0.0, -1.0, 1.0]]).transpose()
        c.init_acc = array([[0.0, 1.0, -1.0]]).transpose()
        c.end_acc = array([[0.0, 100.0, 1.0]]).transpose()
        c.init_jerk = array([[2.0, 4.0, 1.0]]).transpose()
        c.end_jerk = array([[-1.0, 2.0, 7.0]]).transpose()

        c.saveAsText("curve_constraints.txt")
        c.saveAsXML("curve_constraints.xml", "curve_constraints")
        c.saveAsBinary("curve_constraints")

        c_txt = curve_constraints()
        c_txt.loadFromText("curve_constraints.txt")
        self.assertEqual(c, c_txt)

        c_xml = curve_constraints()
        c_xml.loadFromXML("curve_constraints.xml", "curve_constraints")
        self.assertEqual(c, c_xml)

        c_bin = curve_constraints()
        c_bin.loadFromBinary("curve_constraints")
        self.assertEqual(c, c_bin)

        c_pickled = pickle.dumps(c)
        c_from_pickle = pickle.loads(c_pickled)
        self.assertEqual(c_from_pickle, c)


if __name__ == "__main__":
    unittest.main()
