#include<iostream>
#include<cmath>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Rotation_dr( Vector3d &t, Vector3d &tr, double psi, double psi_r )
{
	Matrix3d matrix,I; 
	double sn;
	sn = sin(psi);
	I = Matrix3d::Identity();
	matrix = psi_r*(-sn*I+cos(psi)*cross_vM(t,I))+sn*cross_vM(tr,I);
	return matrix;
}
