#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using namespace Eigen;

void Rotation_dr(const Vector3d & t, const Vector3d & tr, const double psi, const double psi_r,
	               Matrix3d & matrix ) {
	const Matrix3d I = Matrix3d::Identity();
	const double sn = std::sin(psi);

	matrix = psi_r*(-sn*I + std::cos(psi)*cross_vM(t,I))
	       + sn*cross_vM(tr,I);
} // void Rotation_dr(const Vector3d & t, const Vector3d & tr, const double psi, const double psi_r, Matrix3d & matrix ) {
