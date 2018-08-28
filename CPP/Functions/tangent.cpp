#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void tangent(const Vector3d &v, Vector3d & t) {
	t = v/v.norm();
} // void tangent(const Vector3d &v, Vector3d & t) {
