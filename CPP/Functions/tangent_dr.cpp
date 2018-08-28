#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void tangent_dr(const Vector3d & a1, const Vector3d & a1r, const double norm_a1, Vector3d & t) {
	double n = norm_a1;
  t = a1r/n -(a1.dot(a1r))*a1/pow(n,3);
} // void tangent_dr(const Vector3d & a1, const Vector3d & a1r, const double norm_a1, const Vector3d  ) {
