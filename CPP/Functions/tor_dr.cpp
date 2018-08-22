#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
double tor_dr( Matrix3d &LT0T,Matrix3d &LTt,Matrix3d &LT0T1,Matrix3d &LTt1,Matrix3d &LTtr,Matrix3d &LTt1r,Matrix3d &RT,Matrix3d &Rt,Matrix3d &RT1,Matrix3d &Rt1,Matrix3d &Rtr,Matrix3d &Rt1r,Matrix3d &A0,double beta,double alpha )
{
	double scalar;
	Vector3d v1,v2,u1,u2;
	Matrix3d m1,m2,m3,m4;
	
    m1=Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T;
	m2=Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T;
	m3=Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T;
	m4=Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1;
	v1=(m1+m2+m3+m4)*A0.col(beta-1);u1=(Rt*LTt*RT*LT0T)*A0.col(alpha-1);
    
    m1=Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T;
	m2=Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1;
	v2=(m1+m2)*A0.col(beta-1);u2=(Rtr*LTt*RT*LT0T+Rt*LTtr*RT*LT0T)*A0.col(alpha-1);

	scalar = v1.dot(u1) +  v2.dot(u2)   ;
	
	return scalar;
}
