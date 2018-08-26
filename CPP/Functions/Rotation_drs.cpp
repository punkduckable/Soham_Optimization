#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using namespace Eigen;

void Rotation_drs(const Vector3d & t, const Vector3d & tr, const Vector3d & ts, const Vector3d & trs,
	                    const double psi, const double psi_r, const double psi_s, Matrix3d & matrix) {
	// Local variables
	Matrix3d I = Matrix3d::Identity(), CtI = cross_vM(t,I);
	const double sn = std::sin(psi), cs = std::cos(psi), psi_rs = 0;

	matrix = (psi_rs*cs - psi_s*psi_r*sn)*CtI
	       - (psi_rs*sn + psi_s*psi_r*cs)*I
	 			 + (psi_s*cs)*cross_vM(tr,I)
				 + (psi_r*cs)*cross_vM(ts,I)
				 + sn*cross_vM(trs,I);
}
