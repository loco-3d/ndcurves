# Copyright (c) 2020, CNRS
# Authors: Pierre Fernbach <pfernbac@laas.fr>

import unittest

import numpy as np
from numpy import array, array_equal

from ndcurves import constant, constant3


class ConstantCurveTest(unittest.TestCase):
    def test_constructor(self):
        # default constructor
        c = constant()
        self.assertEqual(c.min(), 0.0)
        self.assertEqual(c.max(), 0)
        self.assertEqual(c.dim(), 0)
        self.assertEqual(c.degree(), 0)

        # constructor from a point
        p = array([1, 23.0, 5.0, 9, -5])
        c = constant(p)
        self.assertEqual(c.min(), 0.0)
        self.assertTrue(
            c.max() > 1e100
        )  # convert std::numeric_limits<time_t>::max() to python ?
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 0)

        # constructor with timing
        c = constant(p, 1.0, 5.0)
        self.assertEqual(c.min(), 1.0)
        self.assertEqual(c.max(), 5.0)
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 0)

        with self.assertRaises(ValueError):
            c = constant(p, 2.0, 1.0)

    def test_evaluate(self):
        p = array([1, 23.0, 5.0, 9, -5])
        c = constant(p, 1.0, 3.0)

        self.assertTrue(array_equal(c(1.0), p))
        self.assertTrue(array_equal(c(2.5), p))
        self.assertTrue(array_equal(c(3.0), p))

        with self.assertRaises(ValueError):
            c(0.5)

        with self.assertRaises(ValueError):
            c(4.0)

    def test_derivate(self):
        p = array([1, 23.0, 5.0, 9, -5])
        p0 = np.zeros(5)
        c = constant(p, 1.0, 3.0)

        self.assertTrue(array_equal(c.derivate(1.0, 1), p0))
        self.assertTrue(array_equal(c.derivate(1.0, 2), p0))
        self.assertTrue(array_equal(c.derivate(3.0, 1), p0))
        self.assertTrue(array_equal(c.derivate(1.5, 1), p0))

        with self.assertRaises(ValueError):
            c.derivate(0.5, 1)

        with self.assertRaises(ValueError):
            c.derivate(4.0, 3)

        c_deriv = c.compute_derivate(1)
        self.assertTrue(array_equal(c_deriv(1.5), p0))
        self.assertEqual(c_deriv.min(), c.min())
        self.assertEqual(c_deriv.max(), c.max())
        self.assertEqual(c_deriv.dim(), c.dim())

    def test_comparator(self):
        p = array([1, 23.0, 5.0, 9, -5])
        c1 = constant(p, 1.0, 3.0)
        c2 = constant(p, 1.0, 3.0)
        c3 = c1

        pn = array([1, 23.0, 5.0, 9])
        p2 = array([1, 17.0, 5.0, 9, -5])
        cn1 = constant(pn, 1.0, 3.0)
        cn2 = constant(p, 1.5, 3.0)
        cn3 = constant(p, 1.0, 2.0)
        cn4 = constant(p2, 1.0, 3.0)

        self.assertEqual(c1, c2)
        self.assertEqual(c1, c3)
        self.assertNotEqual(c1, cn1)
        self.assertNotEqual(c1, cn2)
        self.assertNotEqual(c1, cn3)
        self.assertNotEqual(c1, cn4)

        self.assertTrue(c1.isEquivalent(c2))

    def test_serialization(self):
        p = array([1, 23.0, 5.0, 9, -5])
        c = constant(p, 1.0, 3.0)
        c.saveAsText("serialization_curve.txt")
        c.saveAsXML("serialization_curve.xml", "constant")
        c.saveAsBinary("serialization_curve")

        c_txt = constant()
        c_txt.loadFromText("serialization_curve.txt")
        self.assertEqual(c, c_txt)

        c_xml = constant()
        c_xml.loadFromXML("serialization_curve.xml", "constant")
        self.assertEqual(c, c_xml)

        c_bin = constant()
        c_bin.loadFromBinary("serialization_curve")
        self.assertEqual(c, c_bin)


class Constant3CurveTest(unittest.TestCase):
    def test_constructor(self):
        # default constructor
        c = constant()
        self.assertEqual(c.min(), 0.0)
        self.assertEqual(c.max(), 0)
        self.assertEqual(c.dim(), 0)
        self.assertEqual(c.degree(), 0)

        # constructor from a point
        p = array([1, 23.0, 5.0])
        c = constant3(p)
        self.assertEqual(c.min(), 0.0)
        self.assertTrue(
            c.max() > 1e100
        )  # convert std::numeric_limits<time_t>::max() to python ?
        self.assertEqual(c.dim(), 3)
        self.assertEqual(c.degree(), 0)

        # constructor with timing
        c = constant3(p, 1.0, 5.0)
        self.assertEqual(c.min(), 1.0)
        self.assertEqual(c.max(), 5.0)
        self.assertEqual(c.dim(), 3)
        self.assertEqual(c.degree(), 0)

        with self.assertRaises(ValueError):
            c = constant(p, 2.0, 1.0)

    def test_serialization(self):
        p = array([1, 23.0, 5.0])
        c = constant3(p, 0.0, 2.0)
        c.saveAsText("serialization_constant.txt")
        c.saveAsXML("serialization_constant.xml", "constant")
        c.saveAsBinary("serialization_constant")

        c_txt = constant3()
        c_txt.loadFromText("serialization_constant.txt")
        self.assertEqual(c, c_txt)

        c_xml = constant3()
        c_xml.loadFromXML("serialization_constant.xml", "constant")
        self.assertEqual(c, c_xml)

        c_bin = constant3()
        c_bin.loadFromBinary("serialization_constant")
        self.assertEqual(c, c_bin)


if __name__ == "__main__":
    unittest.main()
