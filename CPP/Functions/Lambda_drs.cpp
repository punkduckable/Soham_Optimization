#include<iostream>
#include<cmath>
#include"Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Lambda_drs( Vector3d &T, Vector3d &t, Vector3d &tr, Vector3d &ts, Vector3d &trs )
{
	double dTt; 
	Vector3d CTt,vec_cross1;
	Matrix3d matrix,I,m1,m2,m3,m4,m5,m6,m61,m62; 
	
	dTt=T.dot(t);CTt=T.cross(t);
	I = Matrix3d::Identity(); 
	
    vec_cross1=T.cross(trs);
	
	m1 = T.dot(trs)*I+cross_vM(vec_cross1,I)+(2*T.dot(tr)*T.dot(ts)/pow(1+dTt,3)-T.dot(trs)/pow(1+dTt,2))*CTt*CTt.transpose();
	m2 = -(T.dot(tr)/pow(1+dTt,2))*(CTt*(T.cross(ts)).transpose()+T.cross(ts)*CTt.transpose());
	m3 = -(T.dot(ts)/pow(1+dTt,2))*(CTt*(T.cross(tr)).transpose()+T.cross(tr)*CTt.transpose());
	m4 = (1/(1+dTt))*(CTt*(T.cross(trs)).transpose()+T.cross(ts)*(T.cross(tr)).transpose()+T.cross(tr)*(T.cross(ts)).transpose()+T.cross(trs)*CTt.transpose());
	matrix=m1+m2+m3+m4;
	
	return matrix;
}









