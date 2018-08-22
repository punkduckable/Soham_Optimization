#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_d1r( Vector3d &a1,Vector3d &a11,Vector3d &a1r, Vector3d &a11r, double norm_a1 )
{
	Vector3d t;
	double n = norm_a1;
	t = a11r/n - a1.dot(a1r)*a11/pow(n,3)+3*a1.dot(a11)*a1*a1.dot(a1r)/pow(n,5)-(a1.dot(a11r)+a1r.dot(a11))*a1/pow(n,3)-a1.dot(a11)*a1r/pow(n,3);
	return t;
}
