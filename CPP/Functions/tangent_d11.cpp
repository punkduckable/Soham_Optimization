#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
void tangent_d11(const Vector3d & a1, const Vector3d & a11, const Vector3d & a111, const double norm_a1,
	               Vector3d & t) {
	const double n = norm_a1;
	const double Inv_n = 1/n;
	const double Inv_n_squared = Inv_n*Inv_n;
	const double Inv_n_cubed = Inv_n_squared*Inv_n;
	const double a1_dot_a11 = a1.dot(a11);

	t = a111*Inv_n
	  - a1_dot_a11*a11*Inv_n_cubed
		+ 3*pow(a1_dot_a11,2)*a1*(Inv_n_squared*Inv_n_cubed)
		- (a1.dot(a111) + a11.dot(a11))*a1*Inv_n_cubed
		- a1_dot_a11*a11*Inv_n_cubed;
} // void tangent_d11(const Vector3d & a1, const Vector3d & a11, const Vector3d & a111, const double norm_a1,
