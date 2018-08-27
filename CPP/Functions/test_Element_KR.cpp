#include <iostream>
#include "AllFunctions.cpp"
#include "Eigen/Dense"
#include <time.h>

////////// ROBERT ADDED THE FOLLOWING LINE OF CODE IS FOR TESTING //////////
void Test_Function(void);

using namespace std;
using namespace Eigen;
int main()
{
	// Matrices
	Matrix3d A0;
	MatrixXd P(4,4),Q(4,4);
	P << 0,0.1,0.3,0.7,0,0.3,0.7,0.1,0,0.4,0.9,0.2,0,0.1,0.3,0.5;
	Q << 0.1,0.2,0.7,0.3,0.15,0.5,0.1,0.25,0.3,0.3,0.9,0.55,0,0.1,0.4,0.3;
	A0<<1,0,0,0,1,0,0,0,1;
//	cout << P << endl<< endl;cout << Q << endl<< endl;cout << A0 << endl<< endl;

//	Vectors
	RowVector2d End_pt(2);
	RowVectorXd Ele2Cp(3),knotVector(7),weights(4),Mat(6);
	End_pt<<0,0.5;Ele2Cp <<1,2,3;knotVector <<0,0,0,0.5,1,1,1;weights <<1,1,1,1;Mat<<1,1,1,1,1,1;
//	cout << Ele2Cp << endl<< endl;cout << knotVector << endl<< endl;cout << weights << endl<< endl;cout << Mat << endl<< endl;


	// Scalars
	int ngp, order,nos;
	order=2;ngp=order+8;nos=order+1;
	cout << "ngp: " << ngp << endl << endl;
	MatrixXd K_Element(4*nos,4*nos),F_Element(4*nos,1); // initializing element K(consistent tangent) and R (resiude).

	clock_t timer = clock();

	////////// ROBERT ADDED THIS FOR LOOP FOR TESTING //////////
 	for(unsigned int i = 0; i < 1; i++)
    tie(K_Element,F_Element) = Element_KR(P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,A0,ngp);

	timer = clock() - timer;

	printf("Runtime: %lu ms\n", (1000*timer)/(CLOCKS_PER_SEC));

    cout << "Residue Vector: "<< endl << F_Element << endl << endl;
    cout << "Stiffness Matrix: "<< endl << K_Element << endl << endl;
//    for (int i; i<10; i++)
//    	cout << "Soham" << endl << endl;

	Test_Function();
}

void Test_Function(void) {
	// Needed for testing purposes
	#include <cstdlib>
	#include <ctime>
	#include <cstdio>

	// Seed random number generator using current time.
	srand(time(NULL));

	// Populate random matricies.
	const unsigned long Num_El = 20000000;

	double *s1  = new double[Num_El];
	Vector3d *v1 = new Vector3d[Num_El];
	Matrix3d *m1 = new Matrix3d[Num_El];

	for(unsigned long i = 0; i < Num_El; i++) {
		s1[i] = std::rand();
		v1[i] << std::rand(), std::rand(), std::rand();
		m1[i] = MatrixXd::Random(3,3);
	} // 	for(unsigned long i = 0; i < Num_El; i++) {

	//////////////////////////////////////////////////////////////////////////////
	// Run code many times, time it!

	// Start timer
	clock_t timer = std::clock();

	// Run specified function a lot of times
	for(unsigned long i = 0; i < Num_El - 12; i++)
		s1[i] = tor_dr(m1[i], m1[i+1], m1[i+2], m1[i+3],
			          	 m1[i+4], m1[i+5], m1[i+6], m1[i+7],
									 m1[i+8], m1[i+9], m1[i+10], m1[i+11],
							 		 m1[i+12], 1, 1);

	// Stop timer
	timer = std::clock() - timer;

	// Free dynamically allocated memory
	delete [] s1;
	delete [] v1;
	delete [] m1;

	// Report time
	std::printf("It took %lu ms to run the code %lu times\n", (long)(( 1000*timer )/( (double)CLOCKS_PER_SEC )), Num_El);
}
