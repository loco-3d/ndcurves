# Copyright (c) 2020, CNRS
# Authors: Pierre Fernbach <pfernbac@laas.fr>

import unittest
import NDcurves as curves
from NDcurves import polynomial
import numpy as np
from numpy import array, isclose, array_equal


class MinJerkCurveTest(unittest.TestCase):
    def test_constructors(self):
        # constructor from two points
        init = array([1, 23., 5., 9, -5])
        end = array([1.2, -5.1, 5.6, 1, -2])
        c = polynomial.MinimumJerk(init, end)
        self.assertEqual(c.min(), 0.)
        self.assertEqual(c.max(), 1.)
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 5)
        self.assertTrue(isclose(c(0.), init).all())
        self.assertTrue(isclose(c(1.), end).all())

        # constructor with timing
        c = polynomial.MinimumJerk(init, end, 2., 6.)
        self.assertEqual(c.min(), 2.)
        self.assertEqual(c.max(), 6.)
        self.assertEqual(c.dim(), 5)
        self.assertEqual(c.degree(), 5)
        self.assertTrue(isclose(c(2.), init).all())
        self.assertTrue(isclose(c(6.), end).all())

        end4 = array([1.2, -5.1, 5.6, 1])
        with self.assertRaises(ValueError):
            polynomial.MinimumJerk(init, end4)
        with self.assertRaises(ValueError):
            polynomial.MinimumJerk(init, end, 2., 0.)

    def test_derivate(self):
        init = array([1, 23., 5., 9, -5])
        end = array([1.2, -5.1, 5.6, 1, -2])
        p0 = np.zeros(5)
        c = polynomial.MinimumJerk(init, end, 2., 6.)

        self.assertTrue(isclose(c.derivate(2., 1), p0).all())
        self.assertTrue(isclose(c.derivate(2., 2), p0).all())
        self.assertTrue(isclose(c.derivate(6., 1), p0).all())
        self.assertTrue(isclose(c.derivate(6., 2), p0).all())
        self.assertTrue(isclose(c.derivate(3., 6), p0).all())


if __name__ == '__main__':
    unittest.main()
