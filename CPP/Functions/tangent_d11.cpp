#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_d11( Vector3d &a1,Vector3d &a11, Vector3d &a111, double norm_a1 )
{
	Vector3d t;
	double n = norm_a1;
	t = a111/n - a1.dot(a11)*a11/pow(n,3)+3*pow(a1.dot(a11),2)*a1/pow(n,5)-(a1.dot(a111)+a11.dot(a11))*a1/pow(n,3)-a1.dot(a11)*a11/pow(n,3);
	return t;
}
