#include<iostream>
#include<cmath>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Rotation_drs( Vector3d &t, Vector3d &tr, Vector3d &ts, Vector3d &trs, double psi, double psi_r, double psi_s )
{
	Matrix3d matrix,I,CtI; 
	double sn, cs, psi_rs=0;
	sn = sin(psi);cs = cos(psi);
	I = Matrix3d::Identity();CtI = cross_vM(t,I);
	matrix = psi_rs*(-sn*I+cs*CtI) - psi_s*psi_r*(cs*I+sn*CtI)+psi_s*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(ts,I)+sn*cross_vM(trs,I);
	return matrix;
}
