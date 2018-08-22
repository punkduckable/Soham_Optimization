#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_dr( Vector3d &a1,Vector3d &a1r, double norm_a1 )
{
	Vector3d t;
	double n = norm_a1;
	t = a1r/n -(a1.dot(a1r))*a1/pow(n,3);
	return t;
}
