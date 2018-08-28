#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

void tangent_drs(const Vector3d & a1, const Vector3d & a1r, const Vector3d & a1s, const double norm_a1, Vector3d & t) {
	double n = norm_a1;
	double Inv_n = 1/n;
	double Inv_n_squared = Inv_n*Inv_n;
	double Inv_n_cubed = Inv_n_squared*Inv_n;

	t = - a1r*a1s.dot(a1)*Inv_n_cubed
	  - (a1.dot(a1r)*a1s +a1s.dot(a1r)*a1)*Inv_n_cubed
		+ 3*a1.dot(a1r)*a1*a1.dot(a1s)*(Inv_n_cubed*Inv_n_squared);
} // void tangent_drs(const Vector3d & a1, const Vector3d & a1r, const Vector3d & a1s, const double norm_a1, const Vector3d & t) {
