#ifndef _CLASS_PIECEWISE_CURVE
#define _CLASS_PIECEWISE_CURVE

#include "curve_abc.h"
#include "cubic_hermite_spline.h"
#include "bezier_curve"
#include "polynomial.h"


namespace curves
{
/// \class PiecewiseCurve.
/// 
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
     typename Point= Eigen::Matrix<Numeric, Dim, 1>, 
     typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct piecewise_curve
{



} 

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE