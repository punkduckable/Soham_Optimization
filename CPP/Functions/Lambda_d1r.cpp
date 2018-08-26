#include <iostream>
#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Lambda_d1r(const Vector3d & T, const Vector3d & T1, const Vector3d & t,
	              const Vector3d & t1, const Vector3d & tr, const Vector3d & t1r,
								Matrix3d & matrix) {
	// Local variables
	const double dTt = T.dot(t);
	const Vector3d CTt = T.cross(t);
	const Matrix3d I = Matrix3d::Identity();
	const Vector3d vec_sum1 = (T1.cross(tr) + T.cross(t1r));

	// Intermediate variables. These are used to reduce the quantity and complexity
	// of the computations done to get the result.
	const double Inv_1_plus_dTt = 1/(1 + dTt);
	const double Inv_1_plus_dTt_squared = Inv_1_plus_dTt*Inv_1_plus_dTt;
	const RowVector3d CTt_transpose = CTt.transpose();

	// Calculate components of matrix.
	Matrix3d m1 = (T1.dot(tr) + T.dot(t1r))*I
	   + cross_vM(vec_sum1,I);

	Matrix3d m2 = (2*((T1.dot(t) + T.dot(t1))*T.dot(tr))*Inv_1_plus_dTt*Inv_1_plus_dTt_squared)*CTt*CTt_transpose;

	Matrix3d m3 = (-((T1.dot(tr) + T.dot(t1r))*Inv_1_plus_dTt_squared))*CTt*CTt_transpose;

	Matrix3d m4 = (-((T1.dot(t) + T.dot(t1))*Inv_1_plus_dTt_squared))*(CTt*(T.cross(tr)).transpose()
	                                                       + (T.cross(tr))*CTt_transpose);

	Matrix3d m5 = (-(T.dot(tr)*Inv_1_plus_dTt_squared))*(CTt*(T1.cross(t) + T.cross(t1)).transpose()
	                                         + (T1.cross(t) + T.cross(t1))*CTt_transpose);

	Matrix3d m61 = CTt*(T1.cross(tr) + T.cross(t1r)).transpose()
	    + T.cross(tr)*(T1.cross(t) + T.cross(t1)).transpose();

	Matrix3d m62 = (T1.cross(t) + T.cross(t1))*(T.cross(tr)).transpose()
	    + (T1.cross(tr) + T.cross(t1r))*CTt_transpose;

	Matrix3d m6 = (Inv_1_plus_dTt)*(m61 + m62);

	matrix = m1 + m2 + m3 + m4 + m5 + m6;
} // void Lambda_d1r(const Vector3d & T, const Vector3d & T1, const Vector3d & t, const Vector3d & t1, const Vector3d & tr, const Vector3d & t1r, Matrix3d & matrix) {
