#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

double tor_dr(const Matrix3d & LT0T, const Matrix3d & LTt,   const Matrix3d & LT0T1, const Matrix3d & LTt1,
	            const Matrix3d & LTtr, const Matrix3d & LTt1r, const Matrix3d & RT,    const Matrix3d & Rt,
							const Matrix3d & RT1,  const Matrix3d & Rt1,   const Matrix3d & Rtr,   const Matrix3d & Rt1r,
							const Matrix3d & A0,   const int beta,         const int alpha ) {

	// Calculate v1, u1
	const Matrix3d RT_LT0T = RT*LT0T;
  Matrix3d m1 = (Rt1r*LTt + Rt1*LTtr)*RT_LT0T;

	Matrix3d m2 = (Rtr*LTt1 + Rt*LTt1r)*RT_LT0T;

	const Matrix3d RT1_LT0T = RT1*LT0T;
	const Matrix3d Rtr_LTt_plus_Rt_LTtr = Rtr*LTt + Rt*LTtr;
	Matrix3d m3 = (Rtr_LTt_plus_Rt_LTtr)*RT1_LT0T;

	const Matrix3d RT_LT0T1 = RT*LT0T1;
	Matrix3d m4 = (Rtr_LTt_plus_Rt_LTtr)*RT_LT0T1;

	Vector3d v1 = (m1 + m2 + m3 + m4)*A0.col(beta - 1);
	Vector3d u1 = (Rt*LTt*RT_LT0T)*A0.col(alpha - 1);


	// Calculate v2, u2
  m1 = (Rt1*LTt + Rt*LTt1)*RT_LT0T;

	m2 = Rt*LTt*(RT1_LT0T + RT_LT0T1);

	Vector3d v2 = (m1 + m2)*A0.col(beta-1);
	Vector3d u2 = ((Rtr_LTt_plus_Rt_LTtr)*RT_LT0T)*A0.col(alpha - 1);


	// Calculate result
	return v1.dot(u1) + v2.dot(u2);
} // double tor_dr(const Matrix3d & LT0T, const Matrix3d & LTt,   const Matrix3d & LT0T1, const Matrix3d & LTt1,
