#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

double ben_dr(const Matrix3d & LT0T, const Matrix3d & LTt,   const Matrix3d & LT0T1, const Matrix3d & LTt1,
	            const Matrix3d & LTtr, const Matrix3d & LTt1r, const Matrix3d & RT,    const Matrix3d & Rt,
							const Matrix3d & RT1,  const Matrix3d & Rt1,   const Matrix3d & Rtr,   const Matrix3d & Rt1r,
							const Vector3d & a1,   const Vector3d & a1r,   const Matrix3d & A0,    const int alpha ) {

	// Calculate v1
  Vector3d v1 = ( ( (Rt1r*LTt + Rt1*LTtr + Rtr*LTt1 + Rt*LTt1r)*RT + (Rtr*LTt + Rt*LTtr)*RT1 )*LT0T
		   					+ (Rtr*LTt + Rt*LTtr)*RT*LT0T1
		   					)*A0.col(alpha-1);


	// Calculate v2
	Vector3d v2 = ((Rt1*LTt + Rt*LTt1)*RT*LT0T
			 					+ Rt*LTt*(RT1*LT0T + RT*LT0T1))*A0.col(alpha-1);

//	cout << "v1:" << endl << endl;cout << v1 << endl << endl;
//	cout << "v2:" << endl << endl;cout << v2 << endl << endl;
//
//	cout << "a1:" << endl << endl;cout << a1 << endl << endl;
//	cout << "a1r:" << endl << endl;cout << a1r << endl << endl;
//
//	cout << "A0.col(alpha)" << endl << endl;cout << A0.col(alpha-1) << endl << endl;


	// Calculte reault
	return  v1.dot(a1) + v2.dot(a1r) ;
} // double ben_dr(const Matrix3d & LT0T, const Matrix3d & LTt,   const Matrix3d & LT0T1, const Matrix3d & LTt1,
