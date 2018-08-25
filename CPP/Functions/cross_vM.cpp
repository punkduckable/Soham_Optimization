#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;

Matrix3d cross_vM( const Vector3d & v, const Matrix3d & M ) {
	Matrix3d matrix;

	/* Eigen's matricies are stored in column-major order. To improve performnace,
	we want successive lines of code to modify successive elements of matrix
	(and m) */

	// 0 column
	matrix(0,0) = v(1)*M(2,0) - v(2)*M(1,0);
	matrix(1,0) = v(2)*M(0,0) - v(0)*M(2,0);
	matrix(2,0) = v(0)*M(1,0) - v(1)*M(0,0);

	// 1 colum
  matrix(0,1) = v(1)*M(2,1) - v(2)*M(1,1);
	matrix(1,1) = v(2)*M(0,1) - v(0)*M(2,1);
	matrix(2,1) = v(0)*M(1,1) - v(1)*M(0,1);

	// 2 column
	matrix(0,2) = v(1)*M(2,2) - v(2)*M(1,2);
	matrix(1,2) = v(2)*M(0,2) - v(0)*M(2,2);
	matrix(2,2) = v(0)*M(1,2) - v(1)*M(0,2);

	//	cout << matrix << endl;
	return matrix;
} // Matrix3d cross_vM( Vector3d & v, Matrix3d & M ) {
