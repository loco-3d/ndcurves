/**
* \file Math.h
* \brief Linear algebra and other maths definitions. Based on Eigen 3 or more
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains math definitions used
* used throughout the library.
* Preprocessors definition are used to use eitheir float 
* or double values, and 3 dimensional vectors for
* the Point structure.
*/

#ifndef _SPLINEMATH
#define _SPLINEMATH

#include "Exports.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <deque>
#include <utility> 

namespace spline{

#if (USEFLOAT)
	typedef float Real;
	typedef Eigen::Vector3f Vector3;
	typedef Eigen::Vector2f Vector2;
	typedef Eigen::VectorXf VectorX;
	typedef Eigen::MatrixXf MatrixX;
	typedef Eigen::Matrix4f Matrix4;
	typedef Eigen::Matrix3f Matrix3;
#else
	typedef double Real;
	typedef Eigen::Vector3d Vector3;
	typedef Eigen::Vector2d Vector2;
	typedef Eigen::VectorXd VectorX;
	typedef Eigen::MatrixXd MatrixX;
	typedef Eigen::Matrix4d Matrix4;
	typedef Eigen::Matrix3d Matrix3;
#endif
	
//REF: boulic et al An inverse kinematics architecture enforcing an arbitrary number of strict priority levels
template<typename _Matrix_Type_>
void PseudoInverse(_Matrix_Type_& pinvmat)
{
	Eigen::JacobiSVD<_Matrix_Type_> svd(pinvmat, Eigen::ComputeFullU | Eigen::ComputeFullV);
	VectorX m_sigma = svd.singularValues();

	Real pinvtoler= 1.e-6; // choose your tolerance widely!

	MatrixX m_sigma_inv = MatrixX::Zero(pinvmat.cols(),pinvmat.rows());
	for (long i=0; i<m_sigma.rows(); ++i)
	{
		if (m_sigma(i) > pinvtoler)
			m_sigma_inv(i,i)=1.0/m_sigma(i);
	}
	pinvmat = (svd.matrixV()*m_sigma_inv*svd.matrixU().transpose());
}

} // namespace spline
#endif //_SPLINEMATH
