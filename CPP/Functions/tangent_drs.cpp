#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_drs( Vector3d &a1,Vector3d &a1r, Vector3d &a1s, double norm_a1 )
{
	Vector3d t;
	double n = norm_a1;
	t = - a1r*a1s.dot(a1)/pow(n,3) - (a1.dot(a1r)*a1s +a1s.dot(a1r)*a1)/pow(n,3) + 3*a1.dot(a1r)*a1*a1.dot(a1s)/pow(n,5) ;
	return t;
}
