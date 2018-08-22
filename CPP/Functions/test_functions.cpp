#include <iostream>
#include "AllFunctions.cpp"
#include "Eigen/Dense"
using namespace std;
using namespace Eigen;
int main()
{
    double scalar,alpha,beta;
    double psi, psi_1,psi_r,psi_s,psi_1r,psi_1s,PSI,PSI_1;
    psi = 0.2; psi_1 = 1.5; psi_r = 0.3; psi_s = 0.4; psi_1r = 0.35; psi_1s = 0.45;PSI = 0.9;PSI_1=0.1;
    alpha = 2;beta = 3;

    Vector3d vector,t;
	Vector3d a(2,21,9),a1(1,2,3),ar(11,8,3),as(17,9,27),a11(4,5,6),a1r(45,5,61),a1s(89,51,22),ars(89,51,22);
	Vector3d a111(1,7,8),a11r(13,56,81),a11s(17,35,59),a1rs(19,22,55);
	Vector3d N0(34,56,89),N(43,65,92),A(42,59,19),A0(48,51,79),A1(13,35,52),A01(23,75,12);

	a1=a1/a1.norm();ar=ar/ar.norm();as=as/as.norm();a11=a11/a11.norm();a1r=a1r/a1r.norm();a1s=a1s/a1s.norm();ars=ars/ars.norm();
	a111=a111/a111.norm();a11r=a11r/a11r.norm();a11s=a11s/a11s.norm();a1rs=a1rs/a1rs.norm();
	N0=N0/N0.norm();N=N/N.norm();A=A/A.norm();A0=A0/A0.norm();A1=A1/A1.norm();A01=A01/A01.norm();

	Matrix3d M,matrix;
	Matrix3d LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,Area0;

	M(0,0)=1;M(0,1)=2;M(0,2)=3;M(1,0)=4;M(1,1)=5;M(1,2)=6;M(2,0)=7;M(2,1)=8;M(2,2)=9;
	Area0(0,0)=1;Area0(0,1)=0;Area0(0,2)=0;Area0(1,0)=0;Area0(1,1)=1;Area0(1,2)=0;Area0(2,0)=0;Area0(2,1)=0;Area0(2,2)=1;

//	cout << Area0 << endl;


//	// Check tangent
//	cout <<"Check tangent" << endl;
//	t = tangent(a);
//	cout << t << endl;
//	cout << t.norm()<< endl;
//
//	// Check tangent_d1
//	cout <<"Check tangent_d1" << endl;
//	t = tangent_d1(a1,a11,a1.norm());
//	cout << t << endl;
//
//	cout << a*(a+a).transpose() << endl;
//
//	// Check tangent_dr
//	cout <<"Check tangent_dr" << endl;
//	t = tangent_dr(a1,a1r,a1.norm());
//	cout << t << endl;
//
//	// Check tangent_d11
//	cout <<"Check tangent_d11" << endl;
//	t = tangent_d11(a1,a11,a111,a1.norm());
//	cout << t << endl;
//
//	// Check tangent_d1r
//	cout <<"Check tangent_d11" << endl;
//	t = tangent_d1r(a1,a11,a1r,a11r,a1.norm());
//	cout << t << endl;
//
//	// Check tangent_drs
//	cout <<"Check tangent_drs" << endl;
//	t = tangent_drs(a1,a1r,a1s,a1.norm());
//	cout << t << endl;
//
//	// Check tangent_d1rs
//	cout <<"Check tangent_d1rs" << endl;
//	t = tangent_d1rs(a1,a11,a1r,a1s,a11r,a11s,a1.norm());
//	cout << t << endl;
//
//	// Check cross_VM
//	matrix = cross_vM(a,M);
//	cout <<"Check cross_vM" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation
//	matrix = Rotation(a,psi);
//	cout <<"Check Rotation" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation_d1
//	matrix = Rotation_d1(a,a1,psi,psi_1);
//	cout <<"Check Rotation_d1" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation_dr
//	matrix = Rotation_dr(a,ar,psi,psi_r);
//	cout <<"Check Rotation_dr" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation_d1r
//	matrix = Rotation_d1r(a,a1,ar,a1r,psi,psi_1,psi_r,psi_1r);
//	cout <<"Check Rotation_d1r" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation_drs
//	matrix = Rotation_drs(a,ar,as,ars,psi,psi_r,psi_s);
//	cout <<"Check Rotation_drs" << endl;
//	cout << matrix << endl;
//
//	// Check Rotation_d1rs
//	matrix = Rotation_d1rs(a,a1,ar,as,a1r,a1s,ars,a1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s);
//	cout <<"Check Rotation_d1rs" << endl;
//	cout << matrix << endl;
//
//	// Check Lambda
//	matrix = Lambda( N0,N );
//	cout <<"Check Lambda" << endl;
//	cout << matrix << endl;
//
//
//	// Check Lambda_d1
//	matrix = Lambda_d1( A,A1,a,a1 );
//	cout <<"Check Lambda_d1" << endl;
//	cout << matrix << endl;
//
//	// Check Lambda_dr
//	matrix = Lambda_dr( A,a,ar );
//	cout <<"Check Lambda_dr" << endl;
//	cout << matrix << endl;
//
//	// Check Lambda_d1r
//	matrix = Lambda_d1r( A,A1,a,a1,ar,a1r );
//	cout <<"Check Lambda_d1r" << endl;
//	cout << matrix << endl;
//
//	// Check Lambda_drs
//	matrix = Lambda_drs( A,a,ar,as,ars );
//	cout <<"Check Lambda_drs" << endl;
//	cout << matrix << endl;
//
//	// Check Lambda_d1rs
//	matrix = Lambda_d1rs( A,A1,a,a1,ar,as,a1r,a1s,ars,a1rs );
//	cout <<"Check Lambda_d1rs" << endl;
//	cout << matrix << endl;

	// matrices required for calculating following quantities
	LT0T=Lambda(A0,A);LTt=Lambda(A,a);RT=Rotation(A,PSI);Rt=Rotation(a,psi); // Checked

	LT0T1=Lambda_d1(A0,A01,A,A1);LTt1=Lambda_d1(A,A1,a,a1);RT1=Rotation_d1(A,A1,PSI,PSI_1);Rt1=Rotation_d1(a,a1,psi,psi_1); // Checked

	LTtr=Lambda_dr(A,a,ar);Rtr=Rotation_dr(a,ar,psi,psi_r);LTts=Lambda_dr(A,a,as);Rts=Rotation_dr(a,as,psi,psi_s); // Checked

	LTt1r=Lambda_d1r(A,A1,a,a1,ar,a1r);Rt1r=Rotation_d1r(a,a1,ar,a1r,psi,psi_1,psi_r,psi_1r);// Checked

    LTt1s=Lambda_d1r(A,A1,a,a1,as,a1s);Rt1s=Rotation_d1r(a,a1,as,a1s,psi,psi_1,psi_s,psi_1s);// Checked

	LTtrs=Lambda_drs(A,a,ar,as,ars);Rtrs=Rotation_drs(a,ar,as,ars,psi,psi_r,psi_s);// Checked

	LTt1rs=Lambda_d1rs(A,A1,a,a1,ar,as,ars,a1r,a1s,a1rs);Rt1rs=Rotation_d1rs(a,a1,ar,as,a1r,a1s,ars,a1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s);// Checked

//    cout << LT0T << endl<< endl;cout << LTt << endl<< endl;cout << RT << endl<< endl;cout << Rt << endl<< endl;
//    cout << LT0T1 << endl<< endl;cout << LTt1 << endl<< endl;cout << RT1 << endl<< endl;cout << Rt1 << endl<< endl;
//    cout << LTtr << endl<< endl;cout << Rtr << endl<< endl;cout << LTts << endl<< endl;cout << Rts << endl<< endl;
//    cout << LTt1r << endl<< endl;cout << Rt1r << endl<< endl;cout << LTt1s << endl<< endl;cout << Rt1s << endl<< endl;
//    cout << LTtrs << endl<< endl;cout << Rtrs << endl<< endl;cout << LTt1rs << endl<< endl;cout << Rt1rs << endl<< endl;

	// Check ben_dr
	scalar = ben_dr( LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,a1,a1r,Area0,alpha );
	cout <<"Check ben_dr" << endl;
	cout << scalar << endl;

	// Check ben_drs
	scalar = ben_drs( LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,a1,a1r,a1s,Area0,alpha );
	cout <<"Check ben_drs" << endl;
	cout << scalar << endl;

	// Check tor_dr
	scalar = tor_dr( LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,Area0,beta,alpha );
	cout <<"Check tor_dr" << endl;
	cout << scalar << endl;

	// Check tor_drs
	scalar = tor_drs( LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,Area0,beta,alpha);
	cout <<"Check tor_drs" << endl;
	cout << scalar << endl;

	RowVector3d vector1(1,2,3);cout << vector1 ;


}
