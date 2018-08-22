#include<iostream>
#include<cmath>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Rotation_d1rs( Vector3d &t, Vector3d &t1, Vector3d &tr, Vector3d &ts, Vector3d &t1r, Vector3d &t1s, Vector3d &trs, Vector3d &t1rs, double psi, double psi_1, double psi_r, double psi_s, double psi_1r, double psi_1s )
{
	Matrix3d matrix,I,m1,m2,m3,m4; 
	double sn, cs, psi_rs=0, psi_1rs=0;
	sn = sin(psi);cs = cos(psi); 
	I = Matrix3d::Identity();
	
	m1= (-psi_1rs*sn-psi_1r*psi_s*cs-psi_1s*psi_r*cs-psi_1*psi_rs*cs+psi_1*psi_r*psi_s*sn)*I;
	m2= (psi_1rs*cs-psi_1r*psi_s*sn-psi_1s*psi_r*sn-psi_1*psi_rs*sn-psi_1*psi_r*psi_s*cs)*cross_vM(t,I);
	m3= (psi_1s*cs-psi_1*psi_s*sn)*cross_vM(tr,I)+(psi_1r*cs-psi_1*psi_r*sn)*cross_vM(ts,I)+(psi_rs*cs-psi_r*psi_s*sn)*cross_vM(t1,I);
	m4= (psi_1*cs)*cross_vM(trs,I)+(psi_r*cs)*cross_vM(t1s,I)+(psi_s*cs)*cross_vM(t1r,I)+(sn)*cross_vM(t1rs,I);
	matrix = m1+m2+m3+m4;
	return matrix;
}
