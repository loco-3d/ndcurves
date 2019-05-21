#ifndef _CLASS_PIECEWISE_CURVE
#define _CLASS_PIECEWISE_CURVE

#include <type_traits>

#include "curve_abc.h"
#include "cubic_hermite_spline.h"
#include "bezier_curve"
#include "polynomial.h"
#include "curve_conversion.h"


namespace curves
{
/// \class PiecewiseCurve.
/// 
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
     typename Point= Eigen::Matrix<Numeric, Dim, 1>, 
     typename T_Point= std::vector<Point,Eigen::aligned_allocator<Point>,
     typename Polynomial= polynomial <double, double, 3, true, point_t, t_point_t> > >
struct piecewise_curve
{

	typedef Point 	point_t;
	typedef T_Point t_point_t;
	typedef Time 	time_t;
    typedef Numeric	num_t;

    //typedef polynomial  <double, double, 3, true, point_t, t_point_t> polynomial_t;
    typedef std::vector < Polynomial > t_polynomial_t;
    typedef std::vector<Time> t_vector_time_t;

	public:

	piecewise_curve

	private:
	t_polynomial_t polynomial_curves_; // for curves 0/1/2 : [ curve0, curve1, curve2 ]
	t_vector_time_t time_polynomial_curves_; // for curves 0/1/2 : [ (Tmin0,Tmax0),(Tmin1,Tmax1),(Tmin2,Tmax2) ]
	Numeric size_;
} 

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE