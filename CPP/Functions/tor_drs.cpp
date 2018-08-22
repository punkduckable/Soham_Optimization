#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
double tor_drs( Matrix3d &LT0T,Matrix3d &LTt,Matrix3d &LT0T1,Matrix3d &LTt1,Matrix3d &LTtr,Matrix3d &LTts,Matrix3d &LTt1r,Matrix3d &LTt1s,Matrix3d &LTtrs,Matrix3d &LTt1rs,
Matrix3d &RT,Matrix3d &Rt,Matrix3d &RT1,Matrix3d &Rt1,Matrix3d &Rtr,Matrix3d &Rts,Matrix3d &Rt1r,Matrix3d &Rt1s,Matrix3d &Rtrs,Matrix3d &Rt1rs,
Matrix3d &A0,double beta,double alpha )
{
	double scalar;
	Vector3d v1,v2,u1,u2,v3,v4,u3,u4;
	Matrix3d m1,m2,m3,m4;
	
    m1=Rt1rs*LTt*RT*LT0T+Rt1s*LTtr*RT*LT0T+Rt1r*LTts*RT*LT0T+Rt1*LTtrs*RT*LT0T;
	m2=Rtrs*LTt1*RT*LT0T+Rts*LTt1r*RT*LT0T+Rtr*LTt1s*RT*LT0T+Rt*LTt1rs*RT*LT0T;
	m3=Rtrs*LTt*RT1*LT0T+Rts*LTtr*RT1*LT0T+Rtr*LTts*RT1*LT0T+Rt*LTtrs*RT1*LT0T;
	m4=Rtrs*LTt*RT*LT0T1+Rts*LTtr*RT*LT0T1+Rtr*LTts*RT*LT0T1+Rt*LTtrs*RT*LT0T1;
	v1=(m1+m2+m3+m4)*A0.col(beta-1);u1=(Rt*LTt*RT*LT0T)*A0.col(alpha-1);

    m1=Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T;
	m2=Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T;
	m3=Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T;
	m4=Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1;
	v2=(m1+m2+m3+m4)*A0.col(beta-1);u2=(Rts*LTt*RT*LT0T+Rt*LTts*RT*LT0T)*A0.col(alpha-1);
	
	m1=Rt1s*LTt*RT*LT0T+Rts*LTt1*RT*LT0T+Rt1*LTts*RT*LT0T+Rt*LTt1s*RT*LT0T;
	m2=Rts*LTt*RT1*LT0T+Rts*LTt*RT*LT0T1+Rt*LTts*RT1*LT0T+Rt*LTts*RT*LT0T1;
	v3=(m1+m2)*A0.col(beta-1);u3=(Rtr*LTt*RT*LT0T+Rt*LTtr*RT*LT0T)*A0.col(alpha-1);
    
    m1=Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T;
	m2=Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1;
	v4=(m1+m2)*A0.col(beta-1);u4=(Rtrs*LTt*RT*LT0T+Rts*LTtr*RT*LT0T+Rtr*LTts*RT*LT0T+Rt*LTtrs*RT*LT0T)*A0.col(alpha-1);
 

	scalar = v1.dot(u1) +  v2.dot(u2)  + v3.dot(u3) + v4.dot(u4)  ;
	
	return scalar;
}
