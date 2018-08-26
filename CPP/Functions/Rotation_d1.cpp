#include <iostream>
#include <cmath>
#include "Eigen/Dense"

using namespace Eigen;

void Rotation_d1(const Vector3d & t, const Vector3d & t1, const double psi, const double psi_1,
	               Matrix3d & matrix) {
	const Matrix3d I = Matrix3d::Identity();
	const double sn = std::sin(psi);
//	cout << I << endl;
//  cout << cross_vM(t,I) << endl;
//	cout << cross_vM(t1,I) << endl;

	matrix = psi_1*(-sn*I + std::cos(psi)*cross_vM(t,I))
	       + sn*cross_vM(t1,I);
} // void Rotation_d1( Vector3d & t, Vector3d & t1, double psi, double psi_1 , Matrix3d & matrix) {
