#include<iostream>
#include "../Functions/Eigen/Dense"
using namespace std;
using namespace Eigen;
int FindSpan(int n,int p,double u, RowVectorXd U)
{
/* Determine the knot span index */
/* Algorithm A2.1 The Nurbs Book */
/* Input:
	n: m - p - 1 (Read the book for explanation)
	(m = length of knot vector - 1 )
	p: order of bspline.
	u: parameter at which we want to evaluate bspline basis.
	U: knot vector. */
/* Return:
	mid: the knot span index */

	int low,high,mid;

	if (u == U(n+1) )
	{
		return(n); /* Special case */
	}
	low = p; high = n+1; /* Do binary search */
	mid =(low+high)/2;
	while (u < U(mid) || u >= U(mid+1))
	{
		if (u < U(mid))
		{
	   		high = mid;
		}
		else
		{
			low = mid;
		}
			mid = (low+high)/2;
	}
	return(mid);
}
