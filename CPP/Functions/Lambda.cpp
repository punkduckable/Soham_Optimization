#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void Lambda(const Vector3d & N0, const Vector3d & N, Matrix3d & matrix) {
	// local variables
	const double d = N0.dot(N);
	const Vector3d C = N0.cross(N);
	const Matrix3d I = Matrix3d::Identity();

	// return variable
	matrix = d*I
         + cross_vM(C,I)
         + (1/(1+d))*(C*C.transpose());

} // void Lambda( const Vector3d & N0, const Vector3d & N, Matrix3d & matrix) {
