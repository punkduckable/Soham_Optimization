#include<iostream>
#include<cmath>
#include"Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Lambda_d1( Vector3d &T, Vector3d &T1, Vector3d &t, Vector3d &t1 )
{
	double dTt,dTt1,dT1t; 
	Vector3d CTt,CTt1,CT1t,vec_sum1;
	Matrix3d matrix,I; 
	
	dTt=T.dot(t);dTt1=T.dot(t1);dT1t=T1.dot(t);CTt=T.cross(t);CTt1=T.cross(t1);CT1t=T1.cross(t);
	I = Matrix3d::Identity(); 
	
	vec_sum1 = CTt1+CT1t;
	matrix = (dTt1+dT1t)*I+cross_vM(vec_sum1,I)-((dTt1+dT1t)/pow(1+dTt,2))*CTt*CTt.transpose()+(1/(1+dTt))*(CTt*(CT1t+CTt1).transpose()+(CT1t+CTt1)*CTt.transpose())  ;
	
	return matrix;
}
