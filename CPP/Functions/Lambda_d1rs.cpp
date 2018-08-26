#include <iostream>
#include "Eigen/Dense"
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>

using namespace Eigen;

void Lambda_d1rs(const Vector3d & T,   const Vector3d & T1,   const Vector3d & t,   const Vector3d & t1, 
	               const Vector3d & tr,  const Vector3d & ts,   const Vector3d & t1r, const Vector3d & t1s,
								 const Vector3d & trs, const Vector3d & t1rs, Matrix3d & matrix ) {

	// Local variables
	const double dTt = T.dot(t), dT1t = T1.dot(t), dTt1=T.dot(t1), dTts=T.dot(ts), dTtr=T.dot(tr), dTtrs=T.dot(trs);
	const double dTt1s=T.dot(t1s), dTt1r=T.dot(t1r), dTt1rs=T.dot(t1rs), dT1ts=T1.dot(ts), dT1tr=T1.dot(tr), dT1trs=T1.dot(trs);
	Vector3d CTt, CT1t, CTt1, CTts, CTtr, CTtrs, CTt1s, CTt1r, CTt1rs, CT1ts, CT1tr, CT1trs;
	const Matrix3d I = Matrix3d::Identity();
	Matrix3d m[12];

	// Assign vector variables
	CTt=T.cross(t); CT1t=T1.cross(t); CTt1=T.cross(t1); CTts=T.cross(ts); CTtr=T.cross(tr); CTtrs=T.cross(trs);
	CTt1s=T.cross(t1s); CTt1r=T.cross(t1r); CTt1rs=T.cross(t1rs); CT1ts=T1.cross(ts); CT1tr=T1.cross(tr); CT1trs=T1.cross(trs);

	Vector3d vec_sum1 = (CT1trs + CTt1rs);

	//////////////////////////////////////////////////////////////////////////////
	// Precalculation of quantities that are used multiple times.
	// Here we precalculate several quantities that are used mutliple times in
	// evaluating the various m matricies. This reduces the number of times that
	// we need to perform a particular calculation, thereby reducing the amount of
	// work that the computer must do to execute this function.


	/* this quantity is used in the calculation of m2, m3, m4, m7, and m8 (m[1],
	m[2], m[3], m[6], and m[7]). Precalculating this quantity therefore eliminates
	4 additions. */
	const double dT1t_plus_dTt1 = (dT1t + dTt1);

	/* Division is expensive (it's by far the most expensive of the fundamental
	arithmetic operations). Therefore, we want to minimize the number of times
	that we perform division. In the code below, the quantity 'pow(dTt + 1, 2)' is
	calculated 7 times, while 'pow(dTt + 1, 3) is calculated 5 times.
	Precalculating these two quantities therefore eliminates 10 divisions. We
	replace these divisions with 12 multiplications */
	const double Inv_dTt_plus_1 = 1/(1+dTt);
	const double Inv_dTt_plus_1_squared = Inv_dTt_plus_1*Inv_dTt_plus_1;
	const double Inv_dTt_plus_1_cubed = Inv_dTt_plus_1*Inv_dTt_plus_1_squared;


	/* This quantity is used to calculate m9-m12 (m[8]-m[11]), it is calculated a
	total of 4 times. Precaculating this quantity eliminates 3 vector additions
	and three vector tranposes. */
	const RowVector3d CT1t_plus_CTt1_transpose = (CT1t + CTt1).transpose();

	/* The quantity CTt*CTt.transpos() is used in the calculations for m2, m3,
	and m5. Each time this quantity is calculated we need to perform 1
	matrix-matrix multilication and one transpose, this is wasteful. To improve
	runtime, we can calculate this product once and then use this product to
	calculate m2, m3, and m5 (m[1], m[2], and m[4]). This change eliminates 2
	multiplications and two transposes. */
	const Matrix3d CTt_CTtTranspose = CTt*CTt.transpose();


	//////////////////////////////////////////////////////////////////////////////
	// Calculate result

	m[0] = (dT1trs + dTt1rs)*I
	     + cross_vM(vec_sum1,I);

	m[1] = (-6*((dT1t_plus_dTt1*dTtr*dTts)*Inv_dTt_plus_1_squared*Inv_dTt_plus_1_squared))*CTt_CTtTranspose;

	m[2] = ((2*(dT1ts + dTt1s)*dTtr + 2*dT1t_plus_dTt1*dTtrs)*Inv_dTt_plus_1_cubed)*CTt_CTtTranspose;

	m[3] = (2*dT1t_plus_dTt1*dTtr*Inv_dTt_plus_1_cubed)*(CTt*CTts.transpose() + CTts*CTt.transpose());

	m[4] = (2*(dT1tr + dTt1r)*dTts*Inv_dTt_plus_1_cubed - (dT1trs + dTt1rs)*Inv_dTt_plus_1_squared)*CTt_CTtTranspose;

	m[5] = (-((dT1tr + dTt1r)*Inv_dTt_plus_1_squared))*(CTt*CTts.transpose() + CTts*CTt.transpose());

	m[6] = (2*dT1t_plus_dTt1*dTts*Inv_dTt_plus_1_cubed - (dT1ts + dTt1s)*Inv_dTt_plus_1_squared)*(CTt*CTtr.transpose() + CTtr*CTt.transpose());

	m[7] = (-(dT1t_plus_dTt1*Inv_dTt_plus_1_squared))*(CTt*CTtrs.transpose()
	                                    						 + CTts*CTtr.transpose()
																									 + CTtr*CTts.transpose()
																									 + CTtrs*CTt.transpose());

	m[8] = (2*dTtr*dTts*Inv_dTt_plus_1_cubed - dTtrs*Inv_dTt_plus_1_squared)*(CTt*CT1t_plus_CTt1_transpose + (CT1t + CTt1)*CTt.transpose());

	m[9] = (-(dTtr*Inv_dTt_plus_1_squared))*(CTt*(CT1ts + CTt1s).transpose()
			                                   + CTts*CT1t_plus_CTt1_transpose
															 	 				 + (CT1t + CTt1)*CTts.transpose()
															 	 				 + (CT1ts + CTt1s)*CTt.transpose());

	m[10] = (-(dTts*Inv_dTt_plus_1_squared))*(CTt*(CT1tr + CTt1r).transpose()
	                                       	+ CTtr*CT1t_plus_CTt1_transpose
																 			 	 	+ (CT1t + CTt1)*CTtr.transpose()
																 			    + (CT1tr + CTt1r)*CTt.transpose());

	m[11] = Inv_dTt_plus_1*(CTt*(CT1trs + CTt1rs).transpose()
	                      + CTts*(CT1tr + CTt1r).transpose()
										 		+ CTtr*(CT1ts + CTt1s).transpose()
										 		+ CTtrs*CT1t_plus_CTt1_transpose
										 		+ (CT1t + CTt1)*CTtrs.transpose()
										 		+ (CT1ts + CTt1s)*CTtr.transpose()
										 		+ (CT1tr + CTt1r)*CTts.transpose()
										 		+ (CT1trs + CTt1rs)*CTt.transpose());

	matrix = m[0];
	for(unsigned int i = 1; i < 12; i++)
		matrix += m[i];
}
