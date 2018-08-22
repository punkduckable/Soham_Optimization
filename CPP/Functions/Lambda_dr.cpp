#include<iostream>
#include<cmath>
#include "Eigen/Dense"
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

	matrix = dTtr*I+cross_vM(CTtr,I)-(dTtr/pow(1+dTt,2))*CTt*CTt.transpose()+(1/(1+dTt))*(CTt*CTtr.transpose()+CTtr*CTt.transpose());

	return matrix;
}
