#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_d1rs( Vector3d &a1, Vector3d &a11, Vector3d &a1r, Vector3d &a1s, Vector3d &a11r, Vector3d &a11s, double norm_a1 )
{
	Vector3d t, t1, t2, t3, t4, t5, t6 ;
	double n = norm_a1 ;
	t1 = - a11r*a1s.dot(a1)/pow(n,3) - (a1.dot(a1r)*a11s +a1s.dot(a1r)*a11)/pow(n,3) + 3*a1.dot(a1r)*a11*a1.dot(a1s)/pow(n,5) ;
	t2 = 3*(a1s.dot(a1r)*a1.dot(a11)*a1 + a1.dot(a1r)*(a1s.dot(a11)+a1.dot(a11s))*a1)/pow(n,5) ;
	t3 = 3*a1.dot(a1r)*a1.dot(a11)*a1s/pow(n,5)-15*a1.dot(a1r)*a1.dot(a11)*a1*a1.dot(a1s)/pow(n,7) ;
	t4 = -(a1s.dot(a11r)+a11s.dot(a1r))*a1/pow(n,3);
	t5 = -(a1.dot(a11r)+a11.dot(a1r))*a1s/pow(n,3) + 3*(a1.dot(a11r)+a11.dot(a1r))*a1*a1.dot(a1s)/pow(n,5) ;
	t6 = -(a1s.dot(a11)+a1.dot(a11s))*a1r/pow(n,3) + 3*a1.dot(a11)*a1r*a1.dot(a1s)/pow(n,5) ;

//	cout << t1 << endl << t2 << endl <<t3 << endl <<t4 << endl <<t5 << endl << t6 << endl;

	t = t1 + t2 + t3 + t4 + t5 +  t6  ;
	return t;
}
