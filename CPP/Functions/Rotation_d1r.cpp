#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using namespace Eigen;

void Rotation_d1r(const Vector3d & t, const Vector3d & t1, const Vector3d & tr, const Vector3d & t1r,
	                const double psi,   const double psi_1,  const double psi_r,  const double psi_1r,
									Matrix3d & matrix ) {

	const double sn = std::sin(psi), cs = std::cos(psi);
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d CtI = cross_vM(t,I);

	matrix = (-psi_1r*sn - psi_1*psi_r*cs)*I
	       + (psi_1r*cs - psi_1*psi_r*sn)*CtI
				 + psi_1*cs*cross_vM(tr,I)
				 + psi_r*cs*cross_vM(t1,I)
				 + sn*cross_vM(t1r,I);
} // void Rotation_d1r(Vector3d & t, Vector3d & t1, Vector3d & tr, Vector3d & t1r,
