#include<iostream>
#include<cmath>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Rotation_d1( Vector3d &t, Vector3d &t1, double psi, double psi_1 )
{
	Matrix3d matrix,I; 
	double sn;
	sn = sin(psi);
	I = Matrix3d::Identity();
//	cout << I << endl;
//    cout << cross_vM(t,I) << endl;
//	cout << cross_vM(t1,I) << endl;
	matrix = psi_1*(-sn*I + cos(psi)*cross_vM(t,I)) + sn*cross_vM(t1,I);
	return matrix;
}
