#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Matrix3d cross_vM( Vector3d &v, Matrix3d &M )
{
	Matrix3d matrix;
	matrix(0,0)=v(1)*M(2,0)-v(2)*M(1,0);
    matrix(0,1)=v(1)*M(2,1)-v(2)*M(1,1);
	matrix(0,2)=v(1)*M(2,2)-v(2)*M(1,2);
	matrix(1,0)=v(2)*M(0,0)-v(0)*M(2,0);
	matrix(1,1)=v(2)*M(0,1)-v(0)*M(2,1);
	matrix(1,2)=v(2)*M(0,2)-v(0)*M(2,2);
	matrix(2,0)=v(0)*M(1,0)-v(1)*M(0,0);
	matrix(2,1)=v(0)*M(1,1)-v(1)*M(0,1);
	matrix(2,2)=v(0)*M(1,2)-v(1)*M(0,2);
//	cout << matrix << endl;
	return matrix;
}
