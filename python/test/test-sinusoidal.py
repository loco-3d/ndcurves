# Copyright (c) 2020, CNRS
# Authors: Pierre Fernbach <pfernbac@laas.fr>

import unittest

import numpy as np
from numpy import array, isclose

from ndcurves import sinusoidal


class SinusoidalCurveTest(unittest.TestCase):
    def test_constructor(self):
        # default constructor
        c = sinusoidal()
        self.assertEqual(c.min(), 0.0)
        self.assertEqual(c.max(), 0)
        self.assertEqual(c.dim(), 0)
        self.assertEqual(c.degree(), 1)

        # constructor without time
        p0 = array([1, -23.0, 5.0])
        amp = array([-1, 0.6, 2.8])
        T = 1.5
        phi = 0.2
        c = sinusoidal(p0, amp, T, phi)

        self.assertEqual(c.min(), 0.0)
        self.assertTrue(
            c.max() > 1e100
        )  # convert std::numeric_limits<time_t>::max() to python ?
        self.assertEqual(c.dim(), 3)
        self.assertEqual(c.degree(), 1)

        # constructor with timing
        c = sinusoidal(p0, amp, T, phi, 1.0, 5.0)
        self.assertEqual(c.min(), 1.0)
        self.assertEqual(c.max(), 5.0)
        self.assertEqual(c.dim(), 3)
        self.assertEqual(c.degree(), 1)

        # constructor from stationary points
        p_init = array([1.2, -5, 3.6, 5.0])
        p_end = array([-1.2, -3, 2.7, -3.2])
        c = sinusoidal(1.5, p_init, p_end)
        self.assertEqual(c.min(), 0.0)
        self.assertTrue(
            c.max() > 1e100
        )  # convert std::numeric_limits<time_t>::max() to python ?
        self.assertEqual(c.dim(), 4)
        self.assertEqual(c.degree(), 1)

        # constructor from stationary points with timing
        p_init = array([1.2, -5, 3.6, 5.0])
        p_end = array([-1.2, -3, 2.7, -3.2])
        c = sinusoidal(1.5, p_init, p_end, 1.5, 25)
        self.assertEqual(c.min(), 1.5)
        self.assertEqual(
            c.max(), 25.0
        )  # convert std::numeric_limits<time_t>::max() to python ?
        self.assertEqual(c.dim(), 4)
        self.assertEqual(c.degree(), 1)

        with self.assertRaises(ValueError):
            c = sinusoidal(p0, amp, T, phi, 2.0, 1.0)

    def test_evaluate(self):
        p0 = array([1, -23.0, 5.0])
        amp = array([-1, 0.6, 2.8])
        T = 1.5
        phi = 0.0
        c = sinusoidal(p0, amp, T, phi, 0.0, 20.0)

        self.assertTrue(isclose(c(0), p0).all())
        self.assertTrue(isclose(c(T), p0).all())
        self.assertTrue(isclose(c(T / 2.0), p0).all())
        self.assertTrue(isclose(c(9.0 * T), p0).all())
        self.assertTrue(isclose(c(T / 4.0), p0 + amp).all())
        self.assertTrue(isclose(c(3.0 * T / 4.0), p0 - amp).all())

        with self.assertRaises(ValueError):
            c(-0.5)

        with self.assertRaises(ValueError):
            c(24.0)

    def test_derivate(self):
        p0 = array([1, -23.0, 5.0])
        amp = array([-1, 0.6, 2.8])
        T = 1.5
        phi = 0.0
        c = sinusoidal(p0, amp, T, phi, 0.0, 20.0)

        self.assertTrue(isclose(c.derivate(0.0, 1), (amp * 2.0 * np.pi / T)).all())
        self.assertTrue(isclose(c.derivate(T, 1), (amp * 2.0 * np.pi / T)).all())
        self.assertTrue(isclose(c.derivate(T / 2.0, 1), (-amp * 2.0 * np.pi / T)).all())
        self.assertTrue(isclose(c.derivate(T / 4.0, 1), np.zeros(3)).all())
        self.assertTrue(isclose(c.derivate(3.0 * T / 4.0, 1), np.zeros(3)).all())

        self.assertTrue(
            isclose(
                c.derivate(T / 4.0, 2), (-(amp * 2.0 * np.pi / T) * (2.0 * np.pi / T))
            ).all()
        )
        self.assertTrue(
            isclose(
                c.derivate(3.0 * T / 4.0, 2),
                ((amp * 2.0 * np.pi / T) * (2.0 * np.pi / T)),
            ).all()
        )
        self.assertTrue(isclose(c.derivate(0, 2), np.zeros(3)).all())
        self.assertTrue(isclose(c.derivate(T, 2), np.zeros(3)).all())
        self.assertTrue(isclose(c.derivate(2.0 * T, 2), np.zeros(3)).all())

        with self.assertRaises(ValueError):
            c.derivate(-0.5, 1)

        with self.assertRaises(ValueError):
            c.derivate(54.0, 3)

        with self.assertRaises(OverflowError):
            c.derivate(1.0, -1)

        with self.assertRaises(ValueError):
            c.derivate(1.0, 0)

        for i in range(1, 10):
            c_deriv = c.compute_derivate(i)
            self.assertEqual(c_deriv.min(), c.min())
            self.assertEqual(c_deriv.max(), c.max())
            self.assertEqual(c_deriv.dim(), c.dim())
            self.assertTrue(isclose(c_deriv(0.0), c.derivate(0.0, i)).all())
            self.assertTrue(isclose(c_deriv(0.2), c.derivate(0.2, i)).all())
            self.assertTrue(isclose(c_deriv(0.6), c.derivate(0.6, i)).all())
            self.assertTrue(isclose(c_deriv(0.9), c.derivate(0.9, i)).all())
            self.assertTrue(isclose(c_deriv(T), c.derivate(T, i)).all())
            self.assertTrue(isclose(c_deriv(T * 2.0), c.derivate(T * 2.0, i)).all())

    def test_comparator(self):
        p0 = array([1, -23.0, 5.0])
        amp = array([-1, 0.6, 2.8])
        T = 1.5
        phi = 0.0

        c01 = sinusoidal(p0, amp, T, phi)
        c02 = sinusoidal(p0, amp, T, phi)
        self.assertEqual(c01, c02)

        c1 = sinusoidal(p0, amp, T, phi, 1.0, 10.0)
        c2 = sinusoidal(p0, amp, T, phi, 1.0, 10.0)
        c3 = c1
        self.assertEqual(c1, c2)
        self.assertEqual(c1, c3)
        self.assertTrue(c1.isEquivalent(c2))

        p01 = array([1.5, -23.0, 5.0])
        amp1 = array([-1, -0.6, 2.8])
        cn1 = sinusoidal(p01, amp, T, phi, 1.0, 10.0)
        cn2 = sinusoidal(p0, amp1, 1.4, phi, 1.0, 10.0)
        cn3 = sinusoidal(p0, amp, T, 0.2, 1.0, 10.0)
        cn4 = sinusoidal(p0, amp, T, phi, 0.0, 10.0)
        cn5 = sinusoidal(p0, amp, T, phi, 1.0, 5.0)

        self.assertNotEqual(c1, cn1)
        self.assertNotEqual(c1, cn2)
        self.assertNotEqual(c1, cn3)
        self.assertNotEqual(c1, cn4)
        self.assertNotEqual(c1, cn5)

        self.assertTrue(c1.isEquivalent(c2))

    def test_serialization(self):
        p0 = array([1, -23.0, 5.0])
        amp = array([-1, 0.6, 2.8])
        T = 1.5
        phi = 0.2

        c = sinusoidal(p0, amp, T, phi, 2.0, 5.0)
        c.saveAsText("serialization_sinusoidal.txt")
        c.saveAsXML("serialization_sinusoidal.xml", "sinusoidal")
        c.saveAsBinary("serialization_sinusoidal")

        c_txt = sinusoidal()
        c_txt.loadFromText("serialization_sinusoidal.txt")
        self.assertEqual(c, c_txt)

        c_xml = sinusoidal()
        c_xml.loadFromXML("serialization_sinusoidal.xml", "sinusoidal")
        self.assertEqual(c, c_xml)

        c_bin = sinusoidal()
        c_bin.loadFromBinary("serialization_sinusoidal")
        self.assertEqual(c, c_bin)


if __name__ == "__main__":
    unittest.main()
