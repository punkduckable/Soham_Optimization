#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d Lambda( Vector3d &N0, Vector3d &N )
{
	double d;
	Vector3d C;
	Matrix3d matrix,I; 
	
	d = N0.dot(N); C = N0.cross(N);I = Matrix3d::Identity(); 
	
	matrix = d*I+cross_vM(C,I)+(1/(1+d))*(C*C.transpose()) ;
	
	return matrix;
}
