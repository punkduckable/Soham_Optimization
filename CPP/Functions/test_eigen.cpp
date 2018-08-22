#include <iostream>
#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\AllFunctions.cpp>
#include "Eigen/Dense"
#include <tuple>
using namespace std;
using namespace Eigen;
int main()
{
	Matrix3d M,N;
	Vector3d v;
	int order;
	M=Matrix3d::Zero();
//	cout << M << endl << endl;
//	N.col(0) = M.col(0);
//	cout << N << endl << endl;
	
	VectorXd pts,wts;
	order=10;
	tie(pts,wts) = Gauss_quad(order);
	cout<< pts<< endl<< endl;cout<< wts<< endl<< endl;
		
	int order = 2;	 const double* knot(0,0,0,0.5,1,1,1); double t = 0.25; double* N(1,1,1,1);
}
