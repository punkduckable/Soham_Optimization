#include<iostream>
#include "../Functions/Eigen/Dense"
using namespace std;
using namespace Eigen;
MatrixXd NURBS1DBasis2ndDers(double xi,int order,RowVectorXd knotVector,RowVectorXd weights)
{
    int i,span_index,n,m,der_order;
	// Find the span of knotVector in which requested parameter lies.
//	cout << " Finding span " << endl<< endl;
	m = knotVector.cols()-1;n=m-order-1;
	span_index=FindSpan(n,order,xi,knotVector);

	// Find bspline basis fun and derivative values
//	cout << " Finding bspline " << endl<< endl;
	der_order=2;// highest derivative order requested.
	MatrixXd ders(der_order+1,order+1);// Since we have order+1 basis and derivatives from 0th order till der_order.
	ders=DersBasisFun(span_index,xi,order,der_order,knotVector);

	// Calculate nurbs shape functions and derivatives.
//	cout << " Finding NURBS " << endl<< endl;
	double fac;
	MatrixXd DERS(der_order+1,order+1);

	int    uind   = span_index - order;
	double w = 0.0,dwdxi  = 0.0,d2wdxi = 0.0,wi;
	for(i = 0; i <= order; i++)
    {
        wi = weights(uind+i);
        w  = w + ders(0,i)* wi;
        dwdxi = dwdxi + ders(1,i) * wi;
        d2wdxi = d2wdxi + ders(2,i) * wi;
    }


	for(i = 0; i <= order; i++)
    {
        wi        = weights(uind+i);
        fac       = wi/(w*w);
        DERS(0,i) = ders(0,i)*wi/w;
        DERS(1,i) = (ders(1,i)*w - ders(0,i)*dwdxi) * fac;
        DERS(2,i) = wi*(ders(2,i)/w - 2*ders(1,i)*dwdxi/w/w - ders(0,i)*d2wdxi/w/w + 2*ders(0,i)*dwdxi*dwdxi/w/w/w) ;
    }
    return DERS;
}
