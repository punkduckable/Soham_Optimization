#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void tangent_d1(const Vector3d & a1, const Vector3d & a11, const double norm_a1, Vector3d & t ) {
	const double n = norm_a1;
	t = a11/n
	  - (a1.dot(a11))*a1/std::pow(n,3);
} // void tangent_d1(const Vector3d & a1, const Vector3d & a11, const double norm_a1, Vector3d & t ) {
