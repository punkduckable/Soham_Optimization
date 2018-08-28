#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

double tor_drs(const Matrix3d & LT0T,  const Matrix3d & LTt,    const Matrix3d & LT0T1, const Matrix3d & LTt1,
	             const Matrix3d & LTtr,  const Matrix3d & LTts,   const Matrix3d & LTt1r, const Matrix3d & LTt1s,
							 const Matrix3d & LTtrs, const Matrix3d & LTt1rs, const Matrix3d & RT,    const Matrix3d & Rt,
							 const Matrix3d & RT1,   const Matrix3d & Rt1,    const Matrix3d & Rtr,   const Matrix3d & Rts,
							 const Matrix3d & Rt1r,  const Matrix3d & Rt1s,   const Matrix3d & Rtrs,  const Matrix3d & Rt1rs,
							 const Matrix3d & A0,    const int beta,          const int alpha ) {

	// Local variables
	const Matrix3d RT_LT0T = RT*LT0T;
	const Matrix3d RT1_LT0T = RT1*LT0T;
	const Matrix3d RT_LT0T1 = RT*LT0T1;
	const Matrix3d Rtrs_LTt_plus_Rts_LTtr_plus_Rtr_LTts_plus_Rt_LTtrs = Rtrs*LTt + Rts*LTtr + Rtr*LTts + Rt*LTtrs;
	const Matrix3d Rtr_LTt_plus_Rt_LTtr = Rtr*LTt + Rt*LTtr;
	const Matrix3d Rts_LTt = Rts*LTt;
	const Matrix3d Rt_LTts = Rt*LTts;

	// Calculate v1, u1
  Matrix3d m1 = (Rt1rs*LTt
	    + Rt1s*LTtr
		  + Rt1r*LTts
		  + Rt1*LTtrs)*RT_LT0T;

	Matrix3d m2 = (Rtrs*LTt1
	    + Rts*LTt1r
		  + Rtr*LTt1s
		  + Rt*LTt1rs)*RT_LT0T;

	Matrix3d m3 = (Rtrs_LTt_plus_Rts_LTtr_plus_Rtr_LTts_plus_Rt_LTtrs)*RT1_LT0T;

	Matrix3d m4 = (Rtrs_LTt_plus_Rts_LTtr_plus_Rtr_LTts_plus_Rt_LTtrs)*RT_LT0T1;

	Vector3d v1 = (m1 + m2 + m3 + m4)*A0.col(beta - 1);
	Vector3d u1 = (Rt*LTt*RT_LT0T)*A0.col(alpha - 1);


	// Calculate v2, u2
  m1 = (Rt1r*LTt + Rt1*LTtr)*RT_LT0T;

	m2 = (Rtr*LTt1 + Rt*LTt1r)*RT_LT0T;

	m3 = (Rtr_LTt_plus_Rt_LTtr)*RT1_LT0T;

	m4 = (Rtr_LTt_plus_Rt_LTtr)*RT_LT0T1;

	Vector3d v2 = (m1 + m2 + m3 + m4)*A0.col(beta - 1);
	Vector3d u2 = ((Rts_LTt + Rt_LTts)*RT_LT0T)*A0.col(alpha - 1);


	// Calculate v3, u3
	m1 = (Rt1s*LTt
	    + Rts*LTt1
		  + Rt1*LTts
		  + Rt*LTt1s)*RT_LT0T;

	m2 = Rts_LTt*(RT1_LT0T + RT_LT0T1)
		 + Rt_LTts*(RT1_LT0T + RT_LT0T1);

	Vector3d v3 = (m1 + m2)*A0.col(beta - 1);
	Vector3d u3 = ((Rtr_LTt_plus_Rt_LTtr)*RT_LT0T)*A0.col(alpha - 1);


	// Calculate v4, u4
  m1 = (Rt1*LTt + Rt*LTt1)*RT_LT0T;

	m2 = Rt*LTt*(RT1_LT0T + RT_LT0T1);

	Vector3d v4 = (m1 + m2)*A0.col(beta - 1);
	Vector3d u4 = ((Rtrs_LTt_plus_Rts_LTtr_plus_Rtr_LTts_plus_Rt_LTtrs)*RT_LT0T)*A0.col(alpha-1);


	// Calculate result
	return v1.dot(u1) + v2.dot(u2) + v3.dot(u3) + v4.dot(u4);
} // double tor_drs(const Matrix3d & LT0T,  const Matrix3d & LTt,    const Matrix3d & LT0T1, const Matrix3d & LTt1,
