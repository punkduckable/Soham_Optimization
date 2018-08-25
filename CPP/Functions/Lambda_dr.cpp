#include <iostream>
#include <cmath>
#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Lambda_dr(const Vector3d & T, const Vector3d & t, const Vector3d & tr, Matrix3d & matrix) {
	// Local variables
	double dTt = T.dot(t);
	double dTtr = T.dot(tr);
	Vector3d CTt = T.cross(t);
	Vector3d CTtr = T.cross(tr);
	Matrix3d I = Matrix3d::Identity();

	// Calculated result
	matrix = dTtr*I
	       + cross_vM(CTtr,I)
				 - (dTtr/pow(1 + dTt,2))*CTt*CTt.transpose()
				 + (1/(1+dTt))*(CTt*CTtr.transpose() + CTtr*CTt.transpose());
} // void Lambda_dr(const Vector3d & T, const Vector3d & t, const Vector3d & tr, Matrix3d & matrix) {
