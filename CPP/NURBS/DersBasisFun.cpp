#include<iostream>
#include "../Functions/Eigen/Dense"
using namespace std;
using namespace Eigen;
MatrixXd DersBasisFun(int i,double u,int p,int n,RowVectorXd U)
{
	/* Compute nonzero basis functions and their derivatives.*/
	/* First section is A2.2 modified */
	/* to store functions and knot differences. */
	/* Input:
	    i: knot span index (tells us indices of non zero basis functions)
		u: parameter at which we want to evaluate bspline basis and its derivatives.
		p: order of bspline.
		n: order upto which we want to find the derivatives (n <=p )
		U: knot vector.*/
	/* Output: ders: matrix with on of rows = n and no of columns = p+1 */
	int j,j1,j2,r,k,rk,pk,s1,s2;
	double d,saved,temp;
	MatrixXd ndu(p+1,p+1),ders(n+1,p+1),a(100,100);
	VectorXd left(p+1),right(p+1);

	// Calaulate value of basis functions.
	ndu(0,0)=1.0;
	for (j=1; j<=p; j++)
	{
		left(j) = u-U(i+1-j);
		right(j) = U(i+j)-u;
		saved = 0.0;
		for (r=0; r<j; r++)
		{
			ndu(j,r) = right(r+1)+left(j-r) ;/* Lower triangle */
			temp = ndu(r,j-1)/ndu(j,r);
			ndu(r,j) = saved+right(r+1)*temp; /* Upper triangle */
			saved = left(j-r)*temp;
		}
		ndu(j,j) = saved;
	}

	// Calaulate value of derivatives of basis functions.
	for (j=0; j<=p; j++)
	{ /* Load the basis functions */
		ders (0,j) = ndu (j,p) ;
	}
//	cout << " ders(0,:): "<< ders.row(0)  << endl<< endl;

	/* This section computes the derivatives (Eq. [2.9]) */
	for (r=0; r<=p; r++) /* Loop over function index */
	{
		s1=0; s2=1; /* Alternate rows in array a */
		a(0,0) = 1.0;
		/* Loop to compute kth derivative */
		for (k=1; k<=n; k++)
		{
			d = 0.0;
			rk = r-k; pk = p-k;
			if (r >= k)
			{
				a(s2,0) = a(s1,0)/ndu(pk+1,rk);
				d = a(s2,0)*ndu(rk,pk);
			}
			if (rk >= -1)
			{
				j1 = 1;
			}
			else
			{
				j1 = -rk;
			}
			if (r-1 <= pk)
			{
				j2 = k-1;
			}
			else
			{
				j2 = p-r;
			}
			for (j=j1;j<=j2; j++)
			{
				a(s2,j) = (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j);
				d = d + a(s2,j)*ndu(rk+j,pk);
			}
			if (r <= pk)
			{
				a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);
				d = d + a(s2,k)*ndu(r,pk);
			}
			ders (k,r) = d;
		    j=s1; s1=s2; s2=j; /* Switch rows */
		}
	}
	/* Multiply through by the correct factors */
	/* (Eq. [2.9]) */
	r = p;
	for (k=1; k<=n; k++)
	{
		for (j=0; j<=p; j++)
		{
			ders(k,j) = ders(k,j)*r;
		}
      r = r*(p-k);
	}
	return ders;
}
