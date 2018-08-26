#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

double ben_dr(Matrix3d & LT0T, Matrix3d & LTt,   Matrix3d & LT0T1, Matrix3d & LTt1,
	            Matrix3d & LTtr, Matrix3d & LTt1r, Matrix3d & RT,    Matrix3d & Rt,
							Matrix3d & RT1,  Matrix3d & Rt1,   Matrix3d & Rtr,   Matrix3d & Rt1r,
							Vector3d & a1,   Vector3d & a1r,   Matrix3d & A0,    double alpha ) {

	double scalar;
	Vector3d v1,v2;

  v1 = (Rt1r*LTt*RT*LT0T
		 	+ Rt1*LTtr*RT*LT0T
		 	+ Rtr*LTt1*RT*LT0T
		  + Rt*LTt1r*RT*LT0T
		  + Rtr*LTt*RT1*LT0T
		  + Rt*LTtr*RT1*LT0T
		  + Rtr*LTt*RT*LT0T1
		  + Rt*LTtr*RT*LT0T1)*A0.col(alpha-1);

	v2 = (Rt1*LTt*RT*LT0T
		  + Rt*LTt1*RT*LT0T
			+ Rt*LTt*RT1*LT0T
			+ Rt*LTt*RT*LT0T1)*A0.col(alpha-1);

//	cout << "v1:" << endl << endl;cout << v1 << endl << endl;
//	cout << "v2:" << endl << endl;cout << v2 << endl << endl;
//
//	cout << "a1:" << endl << endl;cout << a1 << endl << endl;
//	cout << "a1r:" << endl << endl;cout << a1r << endl << endl;
//
//	cout << "A0.col(alpha)" << endl << endl;cout << A0.col(alpha-1) << endl << endl;


	scalar = v1.dot(a1) + v2.dot(a1r) ;

	return scalar;
}
