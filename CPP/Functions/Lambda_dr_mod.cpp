#include<iostream>
#include<cmath>
#include"Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Lambda_dr( Vector3d &T, Vector3d &t, Vector3d &tr )
{
	double dTt,dTtr;
	Vector3d CTt,CTtr;
	Matrix3d matrix,I;

	dTt=T.dot(t);dTtr=T.dot(tr);CTt=T.cross(t);CTtr=T.cross(tr);
	I = Matrix3d::Identity();

	double dTt_Plus_1 = dTt + 1;
	double Inv_dTt_Plus_1 = 1./(dTt_Plus_1);
	Matrix3d CTtr_CTtTranspose = CTtr*CTt.transpose();

	matrix = dTtr*I
	       + cross_vM(CTtr,I)
	       - (dTtr*pow(Inv_dTt_Plus_1,2))*CTt*CTt.transpose()
				 + Inv_dTt_Plus_1*(CTtr_CTtTranspose + CTtr_CTtTranspose.transpose());

	return matrix;
}
