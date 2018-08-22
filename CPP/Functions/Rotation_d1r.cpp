#include<iostream>
#include<cmath>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Rotation_d1r( Vector3d &t, Vector3d &t1, Vector3d &tr, Vector3d &t1r, double psi, double psi_1, double psi_r, double psi_1r )
{
	Matrix3d matrix,I,CtI; 
	double sn, cs;
	sn = sin(psi);cs = cos(psi); 
	I = Matrix3d::Identity();CtI = cross_vM(t,I);
	matrix = psi_1r*(-sn*I+cs*CtI)-psi_1*psi_r*(cs*I+sn*CtI)+psi_1*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(t1,I)+sn*cross_vM(t1r,I);
	return matrix;
}
