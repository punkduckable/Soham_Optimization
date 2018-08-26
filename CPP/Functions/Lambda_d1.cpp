#include <iostream>
#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Lambda_d1(const Vector3d & T, const Vector3d & T1, const Vector3d & t, const Vector3d & t1,
	             Matrix3d & matrix) {
	const double dTt = T.dot(t), dTt1 = T.dot(t1), dT1t = T1.dot(t);
	Vector3d CTt, CTt1, CT1t, vec_sum1;
	const Matrix3d I = Matrix3d::Identity();

	CTt = T.cross(t); CTt1 = T.cross(t1); CT1t = T1.cross(t);
	vec_sum1 = CT1t + CTt1;

	matrix = (dTt1 + dT1t)*I
	       + cross_vM((vec_sum1),I)
				 - ((dTt1 + dT1t)/pow(1+dTt,2))*CTt*CTt.transpose()
				 + (1/(1 + dTt))*(CTt*(CT1t + CTt1).transpose()
				               + (CT1t + CTt1)*CTt.transpose());

} // void Lambda_d1(const Vector3d & T, const Vector3d & T1, const Vector3d & t, const Vector3d & t1, Matrix3d & matrix){
