#ifndef _CLASS_CURVE_CONVERSION
#define _CLASS_CURVE_CONVERSION

#include "curve_abc.h"
#include "bernstein.h"
#include "curve_constraint.h"

#include "MathDefs.h"

#include <vector>
#include <stdexcept>

#include <iostream>

namespace curves
{

/// \brief Converts a Bezier curve to a polynomial.
/// \param curve   : the Bezier curve defined between [0,1] to convert.
/// \return the equivalent polynomial.
template<typename Bezier, typename Polynomial>
Polynomial polynom_from_bezier(const Bezier& curve)
{
    typedef typename Polynomial::t_point_t    t_point_t;
    typedef typename Polynomial::num_t    num_t;
    assert (curve.min() == 0.);
    assert (curve.max() == 1.);
    t_point_t coefficients;
    Bezier current (curve);
    coefficients.push_back(curve(0.));
    num_t fact = 1;
    for(std::size_t i = 1; i<= curve.degree_; ++i)
    {
        current = current.compute_derivate(1);
        fact *= (num_t)i;
        coefficients.push_back(current(0.)/fact);
    }
    return Polynomial(coefficients,curve.min(),curve.max());
}

template<typename Hermite, typename Polynomial>
Polynomial polynom_from_hermite(const Hermite& curve)
{
    typedef typename Polynomial::t_point_t    t_point_t;
    typedef typename Polynomial::num_t    num_t;
    assert (curve.min() == 0.);
    assert (curve.max() == 1.);
    t_point_t coefficients;
    Hermite current (curve);
    coefficients.push_back(curve(0.));
    num_t fact = 1;
    for(std::size_t i = 1; i<= 3; ++i)
    {
        fact *= (num_t)i;
        coefficients.push_back(current.derivate(0.,i)/fact);
    }
    return Polynomial(coefficients,curve.min(),curve.max());
}

/// \brief Converts a Cubic Hermite curve to a cubic bezier.
/// \param curve   : the cubic hermite curve defined between [0,1] to convert.
/// \return the equivalent cubic bezier curve.
template<typename Hermite, typename Bezier>
Bezier bezier_from_hermite(const Hermite& curve)
{
	typedef typename Hermite::pair_point_tangent_t pair_point_tangent_t;
	typedef typename Bezier::point_t point_t;
	typedef typename Bezier::t_point_t t_point_t;
    typedef typename Bezier::num_t num_t;
    assert (curve.min() == 0.);
    assert (curve.max() == 1.);

    Hermite current (curve);
    assert(current.control_points_.size() >= 2);
    
    pair_point_tangent_t pair0 = current.control_points_.at(0);
    pair_point_tangent_t pair1 = current.control_points_.at(1);

    // Positions/Velocities of hermite curve
    point_t h_p0 = pair0.first;
    point_t h_m0 = pair0.second;
    point_t h_p1 = pair1.first;
    point_t h_m1 = pair1.second;

    // Convert to bezier control points
    // for t in [0,1] : x'(0)=3(b_p1-b_p0) and x'(1)=3(b_p3-b_p2)
    // so : h_m0=3(b_p1-b_p0) and h_m1=3(b_p3-b_p2)
    // <=> b_p1=(h_m0/3)+b_p0 and b_p2=-(h_m1/3)+b_p3
    point_t b_p0 = h_p0;
    point_t b_p3 = h_p1;
    point_t b_p1 = (h_m0/3)+b_p0;
    point_t b_p2 = -(h_m1/3)+b_p3;

    t_point_t control_points;
    control_points.push_back(b_p0);
    control_points.push_back(b_p1);
    control_points.push_back(b_p2);
    control_points.push_back(b_p3);
    return Bezier(control_points.begin(), control_points.end());
}

} // namespace curve
#endif //_CLASS_CURVE_CONVERSION