#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Lambda_drs(const Vector3d &T, const Vector3d &t, const Vector3d &tr, const Vector3d &ts, const Vector3d &trs, Matrix3d & matrix) {
	// Local variables
	double dTt = T.dot(t);
	Vector3d CTt = T.cross(t);
	Vector3d vec_cross1 = T.cross(trs);
	Matrix3d I = Matrix3d::Identity();
	Matrix3d m1,m2,m3,m4;

	// Intermediate variables (these prevent us from having to calculate the same
	// thing multiple times!)
	double Inv_1_plus_dTt_squared = 1/pow(1 + dTt, 2);

	// Popuate the four m matricies.
	m1 = T.dot(trs)*I
	   + cross_vM(vec_cross1,I)
	   + (2*T.dot(tr)*T.dot(ts)/pow(1 + dTt,3) - T.dot(trs)*Inv_1_plus_dTt_squared)*CTt*CTt.transpose();

	m2 = -(T.dot(tr)*Inv_1_plus_dTt_squared)*(CTt*(T.cross(ts)).transpose()
	                                        + T.cross(ts)*CTt.transpose());

	m3 = -(T.dot(ts)*Inv_1_plus_dTt_squared)*(CTt*(T.cross(tr)).transpose()
	                                				+ T.cross(tr)*CTt.transpose());

	m4 = (1/(1 + dTt))*(CTt*(T.cross(trs)).transpose()
	                  + T.cross(ts)*(T.cross(tr)).transpose()
										+ T.cross(tr)*(T.cross(ts)).transpose()
										+ T.cross(trs)*CTt.transpose());

	matrix = m1;
	matrix += m2;
	matrix += m3;
	matrix += m4;
} // void Lambda_drs(const Vector3d &T, const Vector3d &t, const Vector3d &tr, const Vector3d &ts, const Vector3d &trs, Matrix3d & matrix) {
