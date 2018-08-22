#include<iostream>
#include<cmath>
#include"Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Lambda_d1r( Vector3d &T,Vector3d &T1, Vector3d &t, Vector3d &t1, Vector3d &tr, Vector3d &t1r )
{
	double dTt; 
	Vector3d CTt,vec_sum1;
	Matrix3d matrix,I,m1,m2,m3,m4,m5,m6,m61,m62; 
	
	dTt=T.dot(t);CTt=T.cross(t);
	I = Matrix3d::Identity(); 
	
	vec_sum1=(T1.cross(tr)+T.cross(t1r));
	
	
	m1 = (T1.dot(tr)+T.dot(t1r))*I+cross_vM(vec_sum1,I);
	m2 = 2*(((T1.dot(t)+T.dot(t1))*T.dot(tr))/pow(1+dTt,3))*CTt*CTt.transpose();
	m3 = -((T1.dot(tr)+T.dot(t1r))/pow(1+dTt,2))*CTt*CTt.transpose();
	m4 = -((T1.dot(t)+T.dot(t1))/pow(1+dTt,2))*(CTt*(T.cross(tr)).transpose()+(T.cross(tr))*CTt.transpose());
	m5 = -(T.dot(tr)/pow(1+dTt,2))*(CTt*(T1.cross(t)+T.cross(t1)).transpose()+(T1.cross(t)+T.cross(t1))*CTt.transpose()) ;
	m61 = CTt*(T1.cross(tr)+T.cross(t1r)).transpose()+T.cross(tr)*(T1.cross(t)+T.cross(t1)).transpose();
	m62 = (T1.cross(t)+T.cross(t1))*(T.cross(tr)).transpose()+(T1.cross(tr)+T.cross(t1r))*CTt.transpose();
	m6 = (1/(1+dTt))*(m61+m62) ;
	matrix=m1+m2+m3+m4+m5+m6;
	
	return matrix;
}









