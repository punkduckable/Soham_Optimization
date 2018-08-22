#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent( Vector3d &v )
{
	Vector3d t;
	t = v/v.norm();
	return t;
}
