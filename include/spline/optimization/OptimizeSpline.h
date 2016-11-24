/**
* \file OptimizeSpline.h
* \brief Optimization loop for cubic spline generations
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file uses the mosek library to optimize the waypoints location
* to generate exactCubic spline
*/

#ifndef _CLASS_SPLINEOPTIMIZER
#define _CLASS_SPLINEOPTIMIZER

#include "spline/MathDefs.h" 
#include "spline/exact_cubic.h"

#include "mosek/mosek.h" 
#include <Eigen/SparseCore>

#include <utility>

namespace spline
{
/// \class SplineOptimizer
/// \brief Mosek connection to produce optimized splines
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct SplineOptimizer
{
typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Point point_t;
typedef Time time_t;
typedef Numeric num_t;
typedef exact_cubic<time_t, Numeric, Dim, Safe, Point> exact_cubic_t; 
typedef SplineOptimizer<time_t, Numeric, Dim, Safe, Point> splineOptimizer_t; 

/* Constructors - destructors */
public:
	///\brief Initializes optimizer environment
	SplineOptimizer()
	{
		MSKrescodee  r_ = MSK_makeenv(&env_,NULL);
		assert(r_ == MSK_RES_OK);
	}

	///\brief Destructor
	~SplineOptimizer()
	{
		MSK_deleteenv(&env_);
	}

private:
	SplineOptimizer(const SplineOptimizer&);
	SplineOptimizer& operator=(const SplineOptimizer&);
/* Constructors - destructors */

/*Operations*/
public:
	/// \brief Starts an optimization loop to create curve
	///	\param waypoints : a list comprising at least 2 waypoints in ascending time order
	/// \return An Optimised curve
	template<typename In>
	exact_cubic_t* GenerateOptimizedCurve(In wayPointsBegin, In wayPointsEnd) const;
/*Operations*/

private:
	template<typename In>
	void ComputeHMatrices(In wayPointsBegin, In wayPointsEnd, 
	MatrixX& h1, MatrixX& h2, MatrixX& h3, MatrixX& h4) const;
	
/*Attributes*/
private:
	MSKenv_t env_;
/*Attributes*/
	
private:
	typedef std::pair<time_t, Point> waypoint_t; 
	typedef std::vector<waypoint_t> T_waypoints_t; 
};

template<typename Time, typename Numeric, std::size_t Dim, bool Safe, typename Point>
template<typename In>
inline void SplineOptimizer<Time, Numeric, Dim, Safe, Point>::ComputeHMatrices(In wayPointsBegin, In wayPointsEnd, 
	MatrixX& h1, MatrixX& h2, MatrixX& h3, MatrixX& h4) const
{
	std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
	assert((!Safe) || (size > 1));
	
	In it(wayPointsBegin), next(wayPointsBegin);
	++next;
	Numeric t_previous((*it).first);

	for(std::size_t i(0); next != wayPointsEnd; ++next, ++it, ++i)
	{
		num_t const dTi((*next).first  - (*it).first);
		num_t const dTi_sqr(dTi * dTi);
		// filling matrices values
		h3(i,i)   = -3 / dTi_sqr;
		h3(i,i+1) =  3 / dTi_sqr;
		h4(i,i)   = -2 / dTi;
		h4(i,i+1) = -1 / dTi;
		if( i+2 < size)
		{
			In it2(next); ++ it2;
			num_t const dTi_1((*it2).first - (*next).first);
			num_t const dTi_1sqr(dTi_1 * dTi_1);
			// this can be optimized but let's focus on clarity as long as not needed
			h1(i+1, i)   =  2 / dTi;
			h1(i+1, i+1) =  4 / dTi + 4 / dTi_1;
			h1(i+1, i+2) =  2 / dTi_1;
			h2(i+1, i)   = -6 / dTi_sqr;
			h2(i+1, i+1) = (6 / dTi_1sqr) - (6 / dTi_sqr);
			h2(i+1, i+2) =  6 / dTi_1sqr;
		}
	}
}

template<typename Time, typename Numeric, std::size_t Dim, bool Safe, typename Point>
template<typename In>
inline exact_cubic<Time, Numeric, Dim, Safe, Point>*
	SplineOptimizer<Time, Numeric, Dim, Safe, Point>::GenerateOptimizedCurve(In wayPointsBegin, In wayPointsEnd) const
{
	exact_cubic_t* res = 0;
	int const size((int)std::distance(wayPointsBegin, wayPointsEnd));
	if(Safe && size < 1)
	{
		throw; // TODO
	}
	// refer to the paper to understand all this.
	MatrixX h1 = MatrixX::Zero(size, size);
	MatrixX h2 = MatrixX::Zero(size, size);
	MatrixX h3 = MatrixX::Zero(size, size);
	MatrixX h4 = MatrixX::Zero(size, size);

	// remove this for the time being
	/*MatrixX g1 = MatrixX::Zero(size, size);
	MatrixX g2 = MatrixX::Zero(size, size);
	MatrixX g3 = MatrixX::Zero(size, size);
	MatrixX g4 = MatrixX::Zero(size, size);*/

	ComputeHMatrices(wayPointsBegin, wayPointsEnd, h1, h2, h3, h4);

	// number of Waypoints : T + 1 => T mid points. Dim variables per points, + acceleration + derivations
	// (T * t+ 1 ) * Dim * 3 = nb var 

// NOG const MSKint32t numvar = ( size + size - 1) * 3 * 3;
	const MSKint32t numvar = (size) * Dim * 3;
	/*
	We store the variables in that order to simplifly matrix computation ( see later )
// NOG [ x0  x1 --- xn y0 --- y z0 --- zn x0. --- zn. x0..--- zn.. x0' --- zn..' ] T
	   [ x0  x1 --- xn y0 --- y z0 --- zn x0. --- zn. x0..--- zn..] T
	*/
	
	/*the first constraint is H1x. = H2x => H2x - H1x. = 0
	this will give us size * 3 inequalities constraints
	So this line of A will be writen
	H2 -H1 0 0 0 0
	*/

	int ptOff     = (int) Dim * size; // . offest
	int ptptOff   = (int) Dim * 2 * size; // .. offest
	int prOff     = (int) 3 * ptOff; // ' offest
// NOG int prptOff   = (int) prOff + ptOff; // '. offest
// NOG int prptptOff = (int) prptOff + ptOff; // '.. offest

	MatrixX h2x_h1x = MatrixX::Zero(size * Dim, numvar);
	/**A looks something like that : (n = size)
	[H2(0)  0      0     -H1(0) 0-------------------0]
	[ 0  0  H2(0)  0       0   -H1(0)---------------0]
	[ 0  0  0      H2(0)   0     0     H1(0)--------0]
	...
	[ 0  0  0      0   H2(n)   0    0      0     -H1(n)-0] // row n

	Overall it's fairly easy to fill	
	*/
	for(int i = 0; i < size*Dim; i = i + 3)
	{
		for(int j = 0; j<Dim; ++j)
		{
			int id = i + j;
			h2x_h1x.block(id, j*size    , 1, size)  = h2.row(i%3);
			h2x_h1x.block(id, j*size + ptOff, 1, size) -= h1.row(i%3);
		}
	}

	
	/*the second constraint is x' = G1x + G2x. => G1x + G2x. - x' = 0
	this will give us size * 3 inequalities constraints
	So this line of A will be writen
	H2 -H1 0 0 0 0
	*/MatrixX g1x_g2x = MatrixX::Zero(size * Dim, numvar);
	/**A looks something like that : (n = size)
	[G1(0)  0      0     G2(0) 0-----------------------0 -1 0]
	[ 0  0  G1(0)  0       0   G2(0)-------------------0 -1 0]
	[ 0  0  0      G1(0)   0     0     G2(0)-----------0 -1 0]
	...
	[ 0  0  0      0   G1(n)   0    0      0     G2(n)-0 -1 0] // row n

	Overall it's fairly easy to fill	
	*/
// NOG 
	/*for(int j = 0; j<3; ++j)
	{
		for(int j =0; j<3; ++j)
		{
			int id = i + j;
			g1x_g2x.block(id, j*size	    , 1, size)   = g1.row(i);
			g1x_g2x.block(id, j*size + ptOff, 1, size)   = g2.row(i);
			g1x_g2x.block(id, j*size + prOff, 1, size)  -= MatrixX::Ones(1, size);
		}
	}*/

	/*the third constraint is x.' = G3x + G4x. => G3x + G4x. - x.' = 0
	this will give us size * 3 inequalities constraints
	So this line of A will be writen
	H2 -H1 0 0 0 0
	*/MatrixX g3x_g4x = MatrixX::Zero(size * Dim, numvar);
	/**A looks something like that : (n = size)
	[G3(0)  0      0     G4(0) 0-------------------0 -1 0]
	[ 0  0  G3(0)  0       0   G4(0)---------------0 -1 0]
	[ 0  0  0      G3(0)   0     0     G4(0)--------0 -1 0]
	...
	[ 0  0  0      0   G3(n)   0    0      0     G4(n)-0 -1 0] // row n

	Overall it's fairly easy to fill	
	*/
// NOG 
	/*for(int j = 0; j<3; ++j)
	{
		for(int j =0; j<3; ++j)
		{
			int id = i + j;
			g3x_g4x.block(id, j*size		 , 1, size)   = g3.row(i);
			g3x_g4x.block(id, j*size + ptOff, 1, size)   = g4.row(i);
			g3x_g4x.block(id, j*size + prptOff, 1, size)  -= MatrixX::Ones(1, size);
		}
	}
*/
	/*the fourth constraint is x.. = 1/2(H3x + H4x.) => 1/2(H3x + H4x.) - x.. = 0
	=> H3x + H4x. - 2x.. = 0
	this will give us size * 3 inequalities constraints
	So this line of A will be writen
	H2 -H1 0 0 0 0
	*/
	MatrixX h3x_h4x = MatrixX::Zero(size * Dim, numvar);
	/**A looks something like that : (n = size)
	[H3(0)  0      0     H4(0) 0-------------------0 -2 0]
	[ 0  0  H3(0)  0       0   H4(0)---------------0 -2 0]
	[ 0  0  0      H3(0)   0     0     H4(0)-------0 -2 0]
	...
	[ 0  0  0      0   H3(n)   0    0      0     H4(n)-0 -2 0] // row n

	Overall it's fairly easy to fill	
	*/
	for(int i = 0; i < size*Dim; i = i + Dim)
	{
		for(int j = 0; j<Dim; ++j)
		{
			int id = i + j;
			h3x_h4x.block(id, j*size		   , 1, size) = h3.row(i%3);
			h3x_h4x.block(id, j*size + ptOff  , 1, size) = h4.row(i%3);
			h3x_h4x.block(id, j*size + ptptOff, 1, size) = MatrixX::Ones(1, size) * -2;
		}
	}

	/*the following constraints are easy to understand*/

	/*x0,: = x^0*/
	MatrixX x0_x0 = MatrixX::Zero(Dim, numvar);
	for(int j = 0; j<Dim; ++j)
	{
		x0_x0(0, 0)		   = 1;
		x0_x0(1, size)	   = 1;
		x0_x0(2, size * 2) = 1;
	}

	/*x0.,: = 0*/
	MatrixX x0p_0 = MatrixX::Zero(Dim, numvar);
	for(int j = 0; j<Dim; ++j)
	{
		x0p_0(0, ptOff)	   = 1;
		x0p_0(1, ptOff + size) = 1;
		x0p_0(2, ptOff + size * 2) = 1;
	}
		
	/*xt,: = x^t*/
	MatrixX xt_xt = MatrixX::Zero(Dim, numvar);
	for(int j = 0; j<Dim; ++j)
	{
		xt_xt(0, size - 1)	   = 1;
		xt_xt(1, 2 * size - 1) = 1;
		xt_xt(2, 3* size  - 1) = 1;
	}

	/*xT.,: = 0*/
	MatrixX xtp_0 = MatrixX::Zero(Dim, numvar);
	for(int j = 0; j<Dim; ++j)
	{
		xtp_0(0, ptOff + size - 1)	         = 1;
		xtp_0(1, ptOff + size + size - 1)    = 1;
		xtp_0(2, ptOff + size * 2 + size - 1)= 1;
	}

	//skipping constraints on x and y accelerations for the time being
	// to compute A i'll create an eigen matrix, then i'll convert it to a sparse one and fill those tables

	//total number of constraints
// NOG h2x_h1x (size * Dim) + h3x_h4x (size * Dim ) + g1x_g2x (size * Dim ) + g3x_g4x (size*Dim) 
// NOG + x0_x0 (Dim ) + x0p_0 (Dim) + xt_xt (Dim)  + xtp_0 (Dim) = 4 * Dim * size + 4 * Dim
	// h2x_h1x (size * Dim) + h3x_h4x (size * Dim )
	// + x0_x0 (Dim ) + x0p_0 (Dim) + xt_xt (Dim)  + xtp_0 (Dim) = 2 * Dim * size + 4 * Dim
// NOG const MSKint32t numcon = 12 * size + 12;
	const MSKint32t numcon =  Dim * 2 * size + 4 * Dim; // TODO

	//this gives us the matrix A of size numcon * numvaar
	MatrixX a = MatrixX::Zero(numcon, numvar);
	a.block(0							, 0, size * Dim, numvar) = h2x_h1x;
	a.block(size * Dim					, 0, size * Dim, numvar) = h3x_h4x;
	a.block(size * Dim * 2				, 0,        Dim, numvar) = x0p_0  ;
	a.block(size * Dim * 2 + Dim		, 0,        Dim, numvar) = xtp_0  ;
	a.block(size * Dim * 2 + Dim * 2	, 0,        Dim, numvar) = x0_x0  ;
	a.block(size * Dim * 2 + Dim * 3	, 0,        Dim, numvar) = xt_xt  ;

	//convert to sparse representation
	Eigen::SparseMatrix<Numeric> spA;
	spA = a.sparseView();

	//convert to sparse representation using column
	// http://docs.mosek.com/7.0/capi/Conventions_employed_in_the_API.html#sec-intro-subsubsec-cmo-rmo-matrix

	int nonZeros = spA.nonZeros();
	
	/* Below is the sparse representation of the A
	matrix stored by column. */
	double*		aval  = new double[nonZeros];
	MSKint32t*  asub  = new MSKint32t[nonZeros];
	MSKint32t*  aptrb = new MSKint32t[numvar];
	MSKint32t*  aptre = new MSKint32t[numvar];

	int currentIndex = 0;
	for(int j=0; j<numvar; ++j)
	{
		bool nonZeroAtThisCol = false;
		for(int i=0; i<numcon; ++i)
		{
			if(a(i,j) != 0)
			{
				if(!nonZeroAtThisCol)
				{
					aptrb[j] = currentIndex;
					nonZeroAtThisCol = true;
				}
				aval[currentIndex] = a(i,j);
				asub[currentIndex] = i;
				aptre[j] = currentIndex + 1; //overriding previous value
				++currentIndex;
			}
		}
	}

	/*Q looks like this
	0 0 0 0 0 0  | -> size * 3
	0 0 2 0 0 0  | -> size *3
	0 0 0 0 0 0
	0 0 0 0 2 0
	0 0 0 0 0 0	
	*/
	/* Number of non-zeros in Q.*/
	const MSKint32t  numqz = size * Dim * 2;
	MSKint32t* qsubi = new MSKint32t[numqz]; 
	MSKint32t* qsubj = new MSKint32t[numqz]; 
	double*	   qval  = new double   [numqz];
	for(int id = 0; id < numqz; ++id)
	{
		qsubi[id] = id + ptOff; // we want the x.
		qsubj[id] = id + ptOff;
		 qval[id] = 2;
	}
	
	 /* Bounds on constraints. */
	MSKboundkeye* bkc  = new MSKboundkeye[numcon];
	  
	double*   blc = new double[numcon];
	double*   buc = new double[numcon];

	for(int i = 0; i < numcon - Dim * 2 ; ++i)
	{
		bkc[i] = MSK_BK_FX;
		blc[i] = 0;
		buc[i] = 0;
	}
	for(int i = numcon - Dim * 2; i < numcon - Dim ; ++i) // x0 = x^0
	{
		bkc[i] = MSK_BK_FX;
		blc[i] = wayPointsBegin->second(i - (numcon - Dim * 2) );
		buc[i] = wayPointsBegin->second(i - (numcon - Dim * 2) );
	}
	In last(wayPointsEnd);
	--last;
	for(int i = numcon - 3; i < numcon ; ++i) // xT = x^T
	{
		bkc[i] = MSK_BK_FX;
		blc[i] = last->second(i - (numcon - Dim) );
		buc[i] = last->second(i - (numcon - Dim) );
	}

	///*No Bounds on variables. */
	MSKboundkeye* bkx  = new MSKboundkeye[numvar];
	  
	double*   blx = new double[numvar];
	double*   bux = new double[numvar];

	for(int i = 0; i < numvar; ++i)
	{
		bkx[i] = MSK_BK_FR;
		blx[i] = -MSK_INFINITY;
		bux[i] = +MSK_INFINITY;
	}

	MSKrescodee  r;
	MSKtask_t    task = NULL;
	/* Create the optimization task. */
	r = MSK_maketask(env_,numcon,numvar,&task);

	/* Directs the log task stream to the 'printstr' function. */
	/*if ( r==MSK_RES_OK )
		r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);*/

	/* Append 'numcon' empty constraints.
		The constraints will initially have no bounds. */
	if ( r == MSK_RES_OK )
		r = MSK_appendcons(task,numcon);

	/* Append 'numvar' variables.
		The variables will initially be fixed at zero (x=0). */
	if ( r == MSK_RES_OK )
		r = MSK_appendvars(task,numvar);

	for(int j=0; j<numvar && r == MSK_RES_OK; ++j)
	{
	/* Set the linear term c_j in the objective.*/   
      if(r == MSK_RES_OK) 
        r = MSK_putcj(task,j,0); 

		/* Set the bounds on variable j.
		blx[j] <= x_j <= bux[j] */
		if(r == MSK_RES_OK)
		r = MSK_putvarbound(task,
							j,           /* Index of variable.*/
							bkx[j],      /* Bound key.*/
							blx[j],      /* Numerical value of lower bound.*/
							bux[j]);     /* Numerical value of upper bound.*/

		/* Input column j of A */   
		if(r == MSK_RES_OK)
		r = MSK_putacol(task,
						j,                 /* Variable (column) index.*/
						aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
						asub+aptrb[j],     /* Pointer to row indexes of column j.*/
						aval+aptrb[j]);    /* Pointer to Values of column j.*/
	}

	/* Set the bounds on constraints.
		for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
	for(int i=0; i<numcon && r == MSK_RES_OK; ++i)
		r = MSK_putconbound(task,
							i,           /* Index of constraint.*/
							bkc[i],      /* Bound key.*/
							blc[i],      /* Numerical value of lower bound.*/
							buc[i]);     /* Numerical value of upper bound.*/

	/* Maximize objective function. */
	if (r == MSK_RES_OK)
		r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);


	if ( r==MSK_RES_OK ) 
    {  
    /* Input the Q for the objective. */ 
 
    r = MSK_putqobj(task,numqz,qsubi,qsubj,qval); 
    } 

	if ( r==MSK_RES_OK )
	{
		MSKrescodee trmcode;
    
		/* Run optimizer */
		r = MSK_optimizetrm(task,&trmcode);
		if ( r==MSK_RES_OK )
		{
			double *xx = (double*) calloc(numvar,sizeof(double));
			if ( xx )
			{
				/* Request the interior point solution. */
				MSK_getxx(task,	MSK_SOL_ITR, xx);
				T_waypoints_t nwaypoints;
				In begin(wayPointsBegin);
				for(int i=0; i< size; i = ++i, ++begin)
				{
					point_t target;
					for(int j=0; j< Dim; ++ j)
					{
						target(j) = xx[i + j*size];
					}
					nwaypoints.push_back(std::make_pair(begin->first, target));
				}
				res = new exact_cubic_t(nwaypoints.begin(), nwaypoints.end());
				free(xx);
			}
		}
	}
	/* Delete the task and the associated data. */
	MSK_deletetask(&task);
	return res;
}

} // namespace spline
#endif //_CLASS_SPLINEOPTIMIZER
