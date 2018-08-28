#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void tangent_d1r(const Vector3d & a1, const Vector3d & a11, const Vector3d & a1r, const Vector3d & a11r,
	               const double norm_a1 , Vector3d & t) {
	const double n = norm_a1;
	const double Inv_n = 1/n;
	const double Inv_n_squared = Inv_n*Inv_n;
	const double Inv_n_cubed = Inv_n_squared*Inv_n;

	t = a11r/n
	  - a1.dot(a1r)*a11*(Inv_n_cubed)
	  + 3*a1.dot(a11)*a1*a1.dot(a1r)*(Inv_n_squared*Inv_n_cubed)
		- (a1.dot(a11r)
		+ a1r.dot(a11))*a1*Inv_n_cubed
		- a1.dot(a11)*a1r*Inv_n_cubed;
} // void tangent_d1r(const Vector3d & a1, const Vector3d & a11, const Vector3d & a1r, const Vector3d & a11r,
