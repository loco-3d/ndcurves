# Copyright (c) 2020, CNRS
# Authors: Pierre Fernbach <pfernbac@laas.fr>

import unittest

import numpy as np
from numpy import array, isclose

from ndcurves import polynomial


class MinJerkCurveTest(unittest.TestCase):
    def test_constructors(self):
        # constructor from two points
        init = array([1, 23.0, 5.0, 9, -5])
        end = array([1.2, -5.1, 5.6, 1, -2])
        c = polynomial.MinimumJerk(init, end)
        self.assertEqual(c.min(), 0.0)
        self.assertEqual(c.max(), 1.0)
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 5)
        self.assertTrue(isclose(c(0.0), init).all())
        self.assertTrue(isclose(c(1.0), end).all())

        # constructor with timing
        c = polynomial.MinimumJerk(init, end, 2.0, 6.0)
        self.assertEqual(c.min(), 2.0)
        self.assertEqual(c.max(), 6.0)
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 5)
        self.assertTrue(isclose(c(2.0), init).all())
        self.assertTrue(isclose(c(6.0), end).all())

        end4 = array([1.2, -5.1, 5.6, 1])
        with self.assertRaises(ValueError):
            polynomial.MinimumJerk(init, end4)
        with self.assertRaises(ValueError):
            polynomial.MinimumJerk(init, end, 2.0, 0.0)

    def test_derivate(self):
        init = array([1, 23.0, 5.0, 9, -5])
        end = array([1.2, -5.1, 5.6, 1, -2])
        p0 = np.zeros(5)
        c = polynomial.MinimumJerk(init, end, 2.0, 6.0)

        self.assertTrue(isclose(c.derivate(2.0, 1), p0).all())
        self.assertTrue(isclose(c.derivate(2.0, 2), p0).all())
        self.assertTrue(isclose(c.derivate(6.0, 1), p0).all())
        self.assertTrue(isclose(c.derivate(6.0, 2), p0).all())
        self.assertTrue(isclose(c.derivate(3.0, 6), p0).all())


if __name__ == "__main__":
    unittest.main()
