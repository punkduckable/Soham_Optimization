#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

double tor_drs(const Matrix3d & LT0T,  const Matrix3d & LTt,    const Matrix3d & LT0T1, const Matrix3d & LTt1,
	             const Matrix3d & LTtr,  const Matrix3d & LTts,   const Matrix3d & LTt1r, const Matrix3d & LTt1s,
							 const Matrix3d & LTtrs, const Matrix3d & LTt1rs, const Matrix3d & RT,    const Matrix3d & Rt,
							 const Matrix3d & RT1,   const Matrix3d & Rt1,    const Matrix3d & Rtr,   const Matrix3d & Rts,
							 const Matrix3d & Rt1r,  const Matrix3d & Rt1s,   const Matrix3d & Rtrs,  const Matrix3d & Rt1rs,
							 const Matrix3d & A0,    const double beta,       const double alpha ) {
	double scalar;
	Vector3d v1, v2, u1, u2, v3, v4, u3, u4;
	Matrix3d m1, m2, m3, m4;

	// Calculate v1, u1
  m1 = Rt1rs*LTt*RT*LT0T
	   + Rt1s*LTtr*RT*LT0T
		 + Rt1r*LTts*RT*LT0T
		 + Rt1*LTtrs*RT*LT0T;

	m2 = Rtrs*LTt1*RT*LT0T
	   + Rts*LTt1r*RT*LT0T
		 + Rtr*LTt1s*RT*LT0T
		 + Rt*LTt1rs*RT*LT0T;

	m3 = Rtrs*LTt*RT1*LT0T
	   + Rts*LTtr*RT1*LT0T
		 + Rtr*LTts*RT1*LT0T
		 + Rt*LTtrs*RT1*LT0T;

	m4 = Rtrs*LTt*RT*LT0T1
	   + Rts*LTtr*RT*LT0T1
		 + Rtr*LTts*RT*LT0T1
		 + Rt*LTtrs*RT*LT0T1;

	v1 = (m1 + m2 + m3 + m4)*A0.col(beta - 1);
	u1 = (Rt*LTt*RT*LT0T)*A0.col(alpha - 1);


	// Calculate v2, u2
  m1 = Rt1r*LTt*RT*LT0T
		 + Rt1*LTtr*RT*LT0T;

	m2 = Rtr*LTt1*RT*LT0T
	   + Rt*LTt1r*RT*LT0T;

	m3 = Rtr*LTt*RT1*LT0T
	   + Rt*LTtr*RT1*LT0T;

	m4 = Rtr*LTt*RT*LT0T1
	   + Rt*LTtr*RT*LT0T1;

	v2 = (m1 + m2 + m3 + m4)*A0.col(beta - 1);
	u2 = (Rts*LTt*RT*LT0T
		 + Rt*LTts*RT*LT0T)*A0.col(alpha - 1);


	// Calculate v3, u3
	m1 = Rt1s*LTt*RT*LT0T
	   + Rts*LTt1*RT*LT0T
		 + Rt1*LTts*RT*LT0T
		 + Rt*LTt1s*RT*LT0T;

	m2 = Rts*LTt*RT1*LT0T
		 + Rts*LTt*RT*LT0T1
		 + Rt*LTts*RT1*LT0T
		 + Rt*LTts*RT*LT0T1;

	v3 = (m1 + m2)*A0.col(beta - 1);
	u3 = (Rtr*LTt*RT*LT0T
		 + Rt*LTtr*RT*LT0T)*A0.col(alpha - 1);


	// Calculate v4, u4
  m1 = Rt1*LTt*RT*LT0T
	   + Rt*LTt1*RT*LT0T;

	m2 = Rt*LTt*RT1*LT0T
	   + Rt*LTt*RT*LT0T1;

	v4 = (m1 + m2)*A0.col(beta - 1);
	u4 = (Rtrs*LTt*RT*LT0T
		 + Rts*LTtr*RT*LT0T
		 + Rtr*LTts*RT*LT0T
		 + Rt*LTtrs*RT*LT0T)*A0.col(alpha-1);


	// Calculate result
	scalar = v1.dot(u1) + v2.dot(u2) + v3.dot(u3) + v4.dot(u4);

	return scalar;
}
