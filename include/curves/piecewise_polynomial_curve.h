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
struct piecewise_polynomial_curve
{

	typedef Point 	point_t;
	typedef T_Point t_point_t;
	typedef Time 	time_t;
    typedef Numeric	num_t;

    //typedef polynomial  <double, double, 3, true, point_t, t_point_t> polynomial_t;
    typedef Polynomial polynomial_t;
    typedef std::vector < Polynomial > t_polynomial_t;
    typedef std::vector<Time> t_vector_time_t;

	public:

	piecewise_polynomial_curve(polynomial_t pol)
	{
		size_ = 0;
		T_min_ = pol.min();
		add_polynomial_curve(pol);
	}

	void add_polynomial_curve(polynomial_t pol)
	{
		// Check time continuity : Begin time of pol must be equal to T_max_ of actual piecewise curve.
		if (size_!=0 && pol.min()!=T_max_)
		{
			throw std::out_of_range("Can not add new Polynom to PiecewiseCurve : time discontinuity between T_max_ and pol.min()");
		}
		polynomial_curves_.push_back(pol);
		size_ = polynomial_curves_.size();
		T_max_ = pol.max();
		time_polynomial_curves_.push_back(T_max_);
	}

	private:
	t_polynomial_t polynomial_curves_; // for curves 0/1/2 : [ curve0, curve1, curve2 ]
	t_vector_time_t time_polynomial_curves_; // for curves 0/1/2 : [ Tmin0, Tmax0,Tmax1,Tmax2 ]
	Numeric size_;
	Time T_min_, T_max_;
} 

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE