#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void Lambda( Vector3d &N0, Vector3d &N , Matrix3d & matrix) {
	// local variables
	double d = N0.dot(N);
	Vector3d C = N0.cross(N);
	Matrix3d I = Matrix3d::Identity();

	// return variable
	matrix = d*I
        + cross_vM(C,I)
        + (1/(1+d))*(C*C.transpose());

} // Matrix3d Lambda( Vector3d &N0, Vector3d &N , Matrix3d & matrix) {
