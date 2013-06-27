/**
* \file Exports.h
* \brief Exports definition for the spline dll (windows)
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	#ifdef SPLINE_DLLEXPORT
		#define SPLINE_API __declspec(dllexport)
	#else
		#define SPLINE_API __declspec(dllimport)
	#endif
#else
	#define SPLINE_API 
#endif

