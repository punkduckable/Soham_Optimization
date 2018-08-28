#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

void tangent_d1rs(const Vector3d & a1,   const Vector3d & a11,  const Vector3d & a1r, const Vector3d & a1s,
	                const Vector3d & a11r, const Vector3d & a11s, const double norm_a1, Vector3d & t) {
	// Local variables (these help eliminate division operations)
	double n = norm_a1;
	const double Inv_n = 1/n;
	const double Inv_n_squared = Inv_n*Inv_n;
	const double Inv_n_cubed = Inv_n_squared*Inv_n;
	const double Inv_n_fifth = Inv_n_cubed*Inv_n_squared;

	const double a1_dot_a1r = a1.dot(a1r);
	const double a1s_dot_a1r = a1s.dot(a1r);
	const double a1_dot_a1s = a1.dot(a1s);
	Vector3d t1 = -a11r*a1s.dot(a1)*(Inv_n_cubed)
	            - (a1_dot_a1r*a11s + a1s_dot_a1r*a11)*Inv_n_cubed
		          + 3*a1_dot_a1r*a11*a1_dot_a1s*(Inv_n_cubed*Inv_n_squared);


	const double a1_dot_a11 = a1.dot(a11);
	const double a1s_dot_a11_plus_a1_dot_a11s = a1s.dot(a11) + a1.dot(a11s);
	Vector3d t2 = 3*(a1s_dot_a1r*a1_dot_a11*a1
	            + a1_dot_a1r*(a1s_dot_a11_plus_a1_dot_a11s)*a1)*(Inv_n_fifth);

	Vector3d t3 = 3*a1_dot_a1r*a1_dot_a11*a1s*Inv_n_fifth
	            - 15*a1_dot_a1r*a1_dot_a11*a1*a1_dot_a1s*(Inv_n_fifth*Inv_n_squared);

	Vector3d t4 = -(a1s.dot(a11r)
	            + a11s.dot(a1r))*a1*Inv_n_cubed;

	const double a1_dot_a11r_plus_a11_dot_a1r = a1.dot(a11r) + a11.dot(a1r);
	Vector3d t5 = -(a1_dot_a11r_plus_a11_dot_a1r)*a1s*Inv_n_cubed
	            + 3*(a1_dot_a11r_plus_a11_dot_a1r)*a1*a1_dot_a1s*Inv_n_fifth;

	Vector3d t6 = -(a1s_dot_a11_plus_a1_dot_a11s)*a1r*Inv_n_cubed
	            + 3*a1_dot_a11*a1r*a1_dot_a1s*Inv_n_fifth;

//	cout << t1 << endl << t2 << endl <<t3 << endl <<t4 << endl <<t5 << endl << t6 << endl;

	t = t1 + t2 + t3 + t4 + t5 + t6;
} // void tangent_d1rs(const Vector3d & a1,   const Vector3d & a11,  const Vector3d & a1r, const Vector3d & a1s,
