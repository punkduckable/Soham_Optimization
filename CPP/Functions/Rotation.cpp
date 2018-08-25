#include <iostream>
#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Rotation( const Vector3d & t, const double psi, Matrix3d & matrix) {
	Matrix3d I = Matrix3d::Identity();
//	std::cout << I << endl;
	matrix = I*cos(psi)
	       + sin(psi)*cross_vM(t,I);
} // void Rotation( const Vector3d & t, const double psi, Matrix3d & matrix) {
