#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

double ben_drs(const Matrix3d & LT0T,	 const Matrix3d & LTt,    const Matrix3d & LT0T1,  const Matrix3d & LTt1,
	             const Matrix3d & LTtr,  const Matrix3d & LTts,   const Matrix3d & LTt1r,  const Matrix3d & LTt1s,
							 const Matrix3d & LTtrs, const Matrix3d & LTt1rs, const Matrix3d & RT,     const Matrix3d & Rt,
							 const Matrix3d & RT1, 	 const Matrix3d & Rt1,    const Matrix3d & Rtr,    const Matrix3d & Rts,
							 const Matrix3d & Rt1r,  const Matrix3d & Rt1s,   const Matrix3d & Rtrs,   const Matrix3d & Rt1rs,
							 const Vector3d & a1,    const Vector3d & a1r,    const Vector3d & a1s,    const Matrix3d & A0,
							 const int alpha) {

	// Local variables
	const Matrix3d RT_LT0T = RT*LT0T;
	const Matrix3d RT1_LT0T = RT1*LT0T;
	const Matrix3d RT_LT0T1 = RT*LT0T1;

	// Calculate v1
	Matrix3d m1 = (Rt1rs*LTt
	    				 + Rt1s*LTtr
		  			 	 + Rtrs*LTt1
		           + Rts*LTt1r)*RT_LT0T;

	Matrix3d m2 = (Rtrs*LTt + Rts*LTtr)*RT1_LT0T
		          + (Rtrs*LTt + Rts*LTtr)*RT_LT0T1;

	Matrix3d m3 = (Rt1r*LTts
	             + Rt1*LTtrs
		           + Rtr*LTt1s
		           + Rt*LTt1rs)*RT_LT0T;

	Matrix3d m4 = (Rtr*LTts + Rt*LTtrs)*RT1_LT0T
		          + (Rtr*LTts + Rt*LTtrs)*RT_LT0T1;

	Vector3d v1 = (m1 + m2 + m3 + m4)*A0.col(alpha - 1);


	// Calculate v2
  m1 = (Rt1r*LTt
	    + Rt1*LTtr
		  + Rtr*LTt1
		  + Rt*LTt1r)*RT_LT0T;

	m2 = (Rtr*LTt + Rt*LTtr)*RT1_LT0T
		 + (Rtr*LTt + Rt*LTtr)*RT_LT0T1;

	Vector3d v2 = (m1 + m2)*A0.col(alpha - 1);


	// Calculate v3
	m1 = (Rt1s*LTt
	    + Rt1*LTts
		  + Rts*LTt1
		  + Rt*LTt1s)*RT_LT0T;

	m2 = (Rts*LTt + Rt*LTts)*RT1_LT0T
		 + (Rts*LTt + Rt*LTts)*RT_LT0T1;

	Vector3d v3 = (m1 + m2)*A0.col(alpha - 1);


	// Calculate result
	return v1.dot(a1) + v2.dot(a1s) + v3.dot(a1r);
} // double ben_drs(const Matrix3d & LT0T,	 const Matrix3d & LTt,    const Matrix3d & LT0T1,  const Matrix3d & LTt1,
