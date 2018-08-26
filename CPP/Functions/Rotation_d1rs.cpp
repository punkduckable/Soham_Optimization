#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using namespace Eigen;

void Rotation_d1rs(const Vector3d & t,   const Vector3d & t1,  const Vector3d & tr,  const Vector3d & ts,
	                 const Vector3d & t1r, const Vector3d & t1s, const Vector3d & trs, const Vector3d & t1rs,
								   const double psi,     const double psi_1,   const double psi_r,   const double psi_s,
									 const double psi_1r,  const double psi_1s,  Matrix3d & matrix) {

	Matrix3d m1, m2, m3, m4;
	const Matrix3d I = Matrix3d::Identity();
	const double sn = std::sin(psi), cs = std::cos(psi), psi_rs=0, psi_1rs=0;

	m1 = (-psi_1rs*sn
		  - psi_1r*psi_s*cs
			- psi_1s*psi_r*cs
		  - psi_1*psi_rs*cs
		  + psi_1*psi_r*psi_s*sn)*I;

	m2 = (psi_1rs*cs
		  - psi_1r*psi_s*sn
		  - psi_1s*psi_r*sn
	 	  - psi_1*psi_rs*sn
			- psi_1*psi_r*psi_s*cs)*cross_vM(t,I);

	m3 = (psi_1s*cs - psi_1*psi_s*sn)*cross_vM(tr,I)
	   + (psi_1r*cs - psi_1*psi_r*sn)*cross_vM(ts,I)
		 + (psi_rs*cs - psi_r*psi_s*sn)*cross_vM(t1,I);

	m4 = (psi_1*cs)*cross_vM(trs,I)
	   + (psi_r*cs)*cross_vM(t1s,I)
		 + (psi_s*cs)*cross_vM(t1r,I)
		 + (sn)*cross_vM(t1rs,I);

	matrix = m1;
	matrix += m2;
	matrix += m3;
	matrix += m4;
} // void Rotation_d1rs(Vector3d & t, Vector3d & t1, Vector3d & tr, Vector3d & ts,
