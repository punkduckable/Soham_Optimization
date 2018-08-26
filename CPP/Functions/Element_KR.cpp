#include <iostream>
#include <cmath>
#include"Eigen/Dense"
#include <tuple>
using namespace std;
using namespace Eigen;
 tuple<MatrixXd, MatrixXd>Element_KR(MatrixXd &P, MatrixXd &Q, RowVector2d &End_pt,RowVectorXd &Ele2Cp,RowVectorXd &knotVector,int order,RowVectorXd &weights,RowVectorXd &Mat,Matrix3d &A0,int ngp)
//void Element_KR(MatrixXd &P, MatrixXd &Q, RowVector2d &End_pt,RowVectorXd &Ele2Cp,RowVectorXd &knotVector,int order,RowVectorXd &weights,RowVectorXd &Mat,Matrix3d &A0,int ngp)
{
	// // Element_KR: Gives Consistent Tangent and Residual Vector for an element.

	// //Output

	// F_Element: Vector containing residual vector.
	// Provides residual vector of corresponding dof of all control points in an element.
	// Membrane and Bending contributions considered.

	// // Inputs

	// P: cpt matrix containing x-y-z of all cpts of undeformed configuration. Each column of P has one cpt
	// Q: cpt matrix containing x-y-z of all cpts at some newton iteration. Each column of Q has one cpt
	// End_pt: row vector containing end points (knot vector end pts) of an element.
	// Parent Space : e; Parametric Space : E; Physical Space : x and X.
	// Ele2Cp: a row vector containing index of control point (cpt indices are same as shape fun indices).
	// dof: matrix with degrees of freedom corresponding to the control pts.
	// row: x-y-z components of cpt(no of rows = 3 for 3D) and column = index of cpt (no of columns = total number of cpts).
	// knotVector: Knot vector for a patch to which the element belongs.
	// order: order of NURBS shape functions.
	// weights: wieghts associated with the knot vector.
	// Mat: Material parameters.
	// A0: Refrence triad vectors. (3x3 matrix of those vectors).

	//**************************************************************************************************************//
    int i,dof = 4; // degrees of freedom per cpt. (2 for 2D and 4 for 3D).
    double xi; // parametric point of knotvector corresponding to the gauss points.
	//******** Define intial tangent. ********//
	Vector3d T0,temp_vector;
	temp_vector=A0.col(0);T0=tangent(temp_vector);  // tangent vector to the centreline in refrence triad.
//    cout << T0 <<endl;
    //******** Assign material propeties. ********//
	double Area,E,G,I2,I3,Ip; // material properties and properties of c/s of the beam.
    Area=Mat(0);E=Mat(1);G=Mat(2);I2=Mat(3);I3=Mat(4);Ip=Mat(5);
//	cout << Area <<endl;cout << E<<endl;cout << G<<endl;cout << I2<<endl;cout << I3<<endl;cout << Ip<<endl;
	//******** Form element cpt matrix. ********//
    int nos = order+1; // no of shape functions = no of cpts.
    MatrixXd P_ele(4,nos),Q_ele(4,nos);// P_ele : cpt matrix containing x-y-z of all cpts within element.
	for(i=0; i < nos; i++)
	{
			P_ele.col(i)=P.col(Ele2Cp(i)-1);
	    	Q_ele.col(i)=Q.col(Ele2Cp(i)-1);
	}
//	cout<< P_ele << endl << endl;cout<< Q_ele << endl << endl;
	MatrixXd K_Element,F_Element; // initializing element K(consistent tangent) and R (resiude).
	K_Element = MatrixXd::Zero(4*nos,4*nos);F_Element = MatrixXd::Zero(4*nos,1);
//	cout<< K_Element << endl << endl;cout<< F_Element << endl << endl;
    //******** Get Gauss points and values. ********//
    VectorXd Gauss_e,Gauss_W;
//    cout<< ngp << endl<< endl;
	tie(Gauss_e,Gauss_W) = Gauss_quad(ngp);
//    cout<< Gauss_e << endl<< endl;cout<< Gauss_W << endl<< endl;
	//******** Compute shape fun and derivaties ********//
	MatrixXd N(ngp,nos),dN1(ngp,nos),dN2(ngp,nos);
	double J1=0.5*(End_pt(1)-End_pt(0)); // Jacobian for parent-parameter pair. Note: For 1D no need of determinant.
	int der_order=2; // upto what order derivatives are required.
	MatrixXd ders(der_order+1,order+1);
//	    cout << "Soham" << endl << endl;
	for (i=0; i<ngp; i++)
	{
		xi = 0.5*((End_pt(1)-End_pt(0))*Gauss_e(i)+(End_pt(1)+End_pt(0))); // s is Gauss pt in parametric space.
//		cout<< xi << endl<< endl;
	    ders=NURBS1DBasis2ndDers(xi,order,knotVector,weights);
	    N.row(i)=ders.row(0);dN1.row(i)=ders.row(1);dN2.row(i)=ders.row(2);
//	    cout << "NURBS Basis and derivatives: "<< endl;
//	    cout<< ders << endl<< endl;
	}
//    cout<< N << endl<< endl;cout<< dN1 << endl<< endl;cout<< dN2 << endl<< endl;
//	cout << "Soham" << endl << endl;
	//******** Compute elements of K and R  ********//

	// define required vectors and matrices.
	int j,j1,k,k1,r,s,m1; // variables for for loops and dofs.
	double norm_A1,norm_a1,PSI,PSI_1,psi,psi_1,psi_r,psi_1r,psi_s,psi_1s;
	double Fm,Fb,Ft,Km,Kb,Kt; // Residue and stiffness contributions.
	Vector4d xc,Xc,a1,A1,a11,A11,ar,a1r,a11r,as,a1s,a11s;
	Vector3d T,t,T01,T1,t1,tr,t1r,ts,t1s,trs,t1rs,vec1,vec2,vec3,vec4,vec5,vec6;
	Vector2d Ben,ben,ben_r,Tor,tor,tor_r,ben_s,tor_s,ben_rs,tor_rs;
	Matrix3d LT0T,LTt,RT,Rt,LT0T1,LTt1,RT1,Rt1,LTtr,Rtr,LTt1r,Rt1r,LTts,Rts,LTt1s,Rt1s,LTtrs,Rtrs,LTt1rs,Rt1rs;


	for (i=0;i<ngp;i++) // Gauss point loop
	{
		for (j=0;j<nos;j++) // jth control point
		{
			for (j1=0;j1<dof;j1++) // j1th dof of jth control point.
			{
				r=dof*j+j1; // rth dof (local).
				// Forming required vectors.
        // A1,a1 = tangent vector of the centre line for undeformed and deformed respectively.
        xc=Vector4d::Zero();Xc=xc;a1=xc;A1=xc;a11=xc;A11=xc;ar=xc;a1r=xc;a11r=xc; // set all elements of these to zero.

        for (m1=0;m1<nos;m1++) {
        xc=xc+N(i,m1)*Q_ele.col(m1);Xc=Xc+N(i,m1)*P_ele.col(m1);
        a1=a1+dN1(i,m1)*Q_ele.col(m1);A1=A1+dN1(i,m1)*P_ele.col(m1);
        a11=a11+dN2(i,m1)*Q_ele.col(m1);A11=A11+dN2(i,m1)*P_ele.col(m1);
				}
//        cout << "Forming required vectors: " << endl<< endl;
//        cout << xc << endl<< endl;cout << Xc << endl<< endl;cout << a1 << endl<< endl;cout << A1 << endl<< endl;
//        cout << a11 << endl<< endl;cout << A11 << endl<< endl;

				vec1=A1.segment(0,3);vec2=a1.segment(0,3);
        norm_A1 = vec1.norm();norm_a1 = vec2.norm();
        ar(j1,0)=N(0,j);a1r(j1,0)=dN1(i,j);a11r(j1,0)=dN2(i,j);
        PSI=Xc(3,0);PSI_1=A1(3,0);psi=xc(3,0);psi_1=a1(3,0);
        psi_r=ar(3,0);psi_1r=a1r(3,0);
//        cout << ar << endl<< endl;cout << a1r << endl<< endl;cout << a11r << endl<< endl;

        // Normalized tangent vectors.
        vec1=A1.segment(0,3);T=tangent(vec1);vec1=a1.segment(0,3);t=tangent(vec1);T01=Vector3d::Zero();
        vec1=A1.segment(0,3);vec2=A11.segment(0,3);T1=tangent_d1(vec1,vec2,norm_A1);
        vec1=a1.segment(0,3);vec2=a11.segment(0,3);t1=tangent_d1(vec1,vec2,norm_a1);
        vec1=a1.segment(0,3);vec2=a1r.segment(0,3);tr=tangent_dr(vec1,vec2,norm_a1);
        vec1=a1.segment(0,3);vec2=a11.segment(0,3);vec3=a1r.segment(0,3);vec4=a11r.segment(0,3);t1r=tangent_d1r(vec1,vec2,vec3,vec4,norm_a1);
//        cout << "Normalized tangent vectors: " << endl<< endl;
//        cout << T << endl<< endl;cout << t << endl<< endl;cout << T01 << endl<< endl;
//        cout << T1 << endl<< endl;cout << t1 << endl<< endl;cout << tr << endl<< endl;cout << t1r << endl<< endl;

        // Required Matrices.
        Ben=Vector2d::Zero(); ben=Ben;ben_r=Ben;Tor=Ben; tor=Ben;tor_r=Ben;
        Lambda(T0,T, LT0T);                        Lambda(T,t, LTt);
        Rotation(T,PSI,RT);                        Rotation(t,psi, Rt);
        Lambda_d1(T0,T01,T,T1,LT0T1);              Lambda_d1(T,T1,t,t1, LTt1);
        Rotation_d1(T,T1,PSI,PSI_1,RT1);           Rotation_d1(t,t1,psi,psi_1,Rt1);
        Lambda_dr(T,t,tr, LTtr);                   Rotation_dr(t,tr,psi,psi_r,Rtr);
        Lambda_d1r(T,T1,t,t1,tr,t1r,LTt1r);        Rotation_d1r(t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r,Rt1r);
//        cout << "Required Matrices: " << endl<< endl;
//        cout << " LT0T: " << endl << LT0T << endl << endl;cout << " LTt: " << endl << LTt << endl << endl;
//        cout << " RT: " << endl << RT << endl << endl;cout << " Rt: " << endl << Rt << endl << endl;
//        cout << " LT0T1: " << endl << LT0T1 << endl << endl;cout << " LTt1: " << endl << LTt1 << endl << endl;
//        cout << " RT1: " << endl << RT1 << endl << endl;cout << " Rt1: " << endl << Rt1 << endl << endl;
//        cout << " LTtr: " << endl << LTtr << endl << endl;cout << " Rtr: " << endl << Rtr << endl << endl;
//        cout << " LTt1r: " << endl << LTt1r << endl << endl;cout << " Rt1r: " << endl << Rt1r << endl << endl;

        // Bending terms in Residue
        for(m1 = 0; m1 < 2; m1++) {
          vec1=A1.segment(0,3);Ben(m1)=vec1.dot((RT1*LT0T+RT*LT0T1)*A0.col(m1+1));
          vec1=a1.segment(0,3);ben(m1)=vec1.dot((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0.col(m1+1));
          vec1=a1.segment(0,3);vec2=a1r.segment(0,3);
			    ben_r(m1) = ben_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,vec1,vec2,A0,m1+2); // last argument is alpha/beta (look into the function for more detail)
				}
//        cout << "Bending terms in Residue: " << endl<< endl;
//        cout << Ben<<endl<< endl;cout << ben<<endl<< endl;cout << ben_r<<endl<< endl;

				// Torsion terms in Residue
				vec1=(RT1*LT0T+RT*LT0T1)*A0.col(2);vec2=RT*LT0T*A0.col(1);
				Tor(0)= vec1.dot(vec2);
				vec1=(Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0.col(2);vec2=Rt*LTt*RT*LT0T*A0.col(1);
				tor(0)= vec1.dot(vec2);
				tor_r(0)= tor_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,A0,3,2);
				vec1=(RT1*LT0T+RT*LT0T1)*A0.col(1);vec2=RT*LT0T*A0.col(2);
				Tor(1)= vec1.dot(vec2);
				vec1=(Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0.col(1);vec2=Rt*LTt*RT*LT0T*A0.col(2);
        tor(1)= vec1.dot(vec2);
        tor_r(1)= tor_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,A0,2,3);
//        cout << "Torsion terms in Residue: " << endl<< endl;
//        cout << Tor<<endl<< endl;cout << tor<<endl<< endl;cout << tor_r<<endl<< endl;

        // Membrane Contribution.
        vec1=a1.segment(0,3);vec2=A1.segment(0,3);vec3=a1r.segment(0,3);
        Fm = (Area*E*0.5*(vec1.dot(vec1)-vec2.dot(vec2))*vec1.dot(vec3))/pow(vec2.norm(),3);
//        cout << "Fm" << endl << Fm << endl<< endl;

        // Bending Contribution.
        vec2=A1.segment(0,3);
        Fb = (E/pow(vec2.norm(),3))*(I3*(ben(0)-Ben(0))*ben_r(0)+I2*(ben(1)-Ben(1))*ben_r(1));
//        cout << "Fb" << endl << Fb << endl<< endl;

        // Torsion Contribution.
        vec2=A1.segment(0,3);
        Ft = (0.5*G*Ip/vec2.norm())*((tor(0)-Tor(0))*tor_r(0)+(tor(1)-Tor(1))*tor_r(1));
//        cout << "Ft" << endl << Ft << endl<< endl;

        // Total element residue vector.
        F_Element(r)=F_Element(r)+(Fm+Fb+Ft)*J1*Gauss_W(i);

        for (k=0;k<nos;k++) {
        	for (k1=0;k1<dof;k1++) {
        		s=dof*k + k1;// sth dof (local).

            if (r<=s) {
              // Forming required vectors.
              as = Vector4d::Zero(); a1s = as; a11s = as;
              as(k1,0) = N(0,k); a1s(k1,0) = dN1(i,k); a11s(k1,0) = dN2(i,k); psi_s = as(3,0); psi_1s = a1s(3,0);

              //Normalized tangent vectors.
              vec1 = a1.segment(0,3);   vec2 = a1s.segment(0,3);    ts = tangent_dr(vec1,vec2,norm_a1);
              vec1 = a1.segment(0,3);   vec2 = a11.segment(0,3);    vec3 = a1s.segment(0,3);    vec4 = a11s.segment(0,3);   t1s = tangent_d1r(vec1,vec2,vec3,vec4,norm_a1);
              vec3 = a1r.segment(0,3);  vec4 = a1s.segment(0,3);    trs = tangent_drs(vec1,vec3,vec4,norm_a1);
              vec1 = a1.segment(0,3);   vec2 = a11.segment(0,3);    vec3 = a1r.segment(0,3);    vec4 = a1s.segment(0,3);    vec5 = a11r.segment(0,3);   vec6 = a11s.segment(0,3);
              t1rs = tangent_d1rs(vec1,vec2,vec3,vec4,vec5,vec6,norm_a1);
//              cout << "Normalized tangent vectors: " << endl<< endl;
//              cout << ts << endl<< endl;cout << t1s << endl<< endl;cout << trs << endl<< endl;cout << t1rs << endl<< endl;

              // Required Matrices.
              ben_s = Vector2d::Zero();                              tor_s = ben_s; ben_rs = ben_s; tor_rs = ben_s;
              Lambda_dr(T,t,ts,LTts);                                Rotation_dr(t,ts,psi,psi_s,Rts);
              Lambda_d1r(T,T1,t,t1,ts,t1s,LTt1s);                    Rotation_d1r(t,t1,ts,t1s,psi,psi_1,psi_s,psi_1s, Rt1s);
              Lambda_drs(T,t,tr,ts,trs,LTtrs);                       Rotation_drs(t,tr,ts,trs,psi,psi_r,psi_s,Rtrs);
              Lambda_d1rs(T,T1,t,t1,tr,ts,t1r,t1s,trs,t1rs,LTt1rs);  Rotation_d1rs(t,t1,tr,ts,t1r,t1s,trs,t1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s,Rt1rs);
//              cout << " LTts: " << endl << LTts << endl << endl;cout << " Rts: " << endl << Rts << endl << endl;
//              cout << " LTtrs: " << endl << LTtrs << endl << endl;cout << " Rtrs: " << endl << Rtrs << endl << endl;
//              cout << " LTt1s: " << endl << LTt1s << endl << endl;cout << " Rt1s: " << endl << Rt1s << endl << endl;
//              cout << " LTt1rs: " << endl << LTt1rs << endl << endl;cout << " Rt1rs: " << endl << Rt1rs << endl << endl;

              // Bending terms in Stiffness
              for(m1 = 0; m1 < 2; m1++) {
                vec1 = a1.segment(0,3); vec2 = a1r.segment(0,3); vec3 = a1s.segment(0,3);
				        ben_s(m1) = ben_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,vec1,vec3,A0,m1+2); // last argument is alpha/beta (look into the function for more detail)
							  ben_rs(m1) = ben_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,vec1,vec2,vec3,A0,m1+2);
						  } // for(m1 = 0; m1 < 2; m1++) {
                
//              cout << "Bending terms in Stiffness: " << endl<< endl;
//              cout << ben_s <<endl<< endl;cout <<ben_rs <<endl<< endl;

              // Torsion terms in Stiffness
              tor_s(0) = tor_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,A0,3,2);
              tor_rs(0) = tor_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,A0,3,2);
              tor_s(1) = tor_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,A0,2,3);
              tor_rs(1) = tor_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,A0,2,3);
//              cout << "Torsion terms in Stiffness: " << endl<< endl;
//              cout << tor_s <<endl<< endl;cout << tor_rs <<endl<< endl;

              //Membrane Contribution
              vec1 = a1.segment(0,3); vec2 = a1r.segment(0,3); vec3 = a1s.segment(0,3); vec4 = A1.segment(0,3);
              Km = Area*E*(vec1.dot(vec3)*vec1.dot(vec2)+ 0.5*(vec1.squaredNorm() - vec4.squaredNorm())*vec2.dot(vec3))/pow(vec4.norm(),3);
//              cout << "Km" << endl << Km << endl<< endl;

              //Bending Contribution
              Kb = (E/pow(vec4.norm(),3))*(I3*ben_s(0)*ben_r(0)
                                         + I2*ben_s(1)*ben_r(1)
                                         + I3*(ben(0)-Ben(0))*ben_rs(0)
                                         + I2*(ben(1)-Ben(1))*ben_rs(1));
//              cout << "Kb" << endl <<  Kb << endl<< endl;

              //Torsion Contribution
              Kt = (0.5*G*Ip/vec4.norm())*(tor_r(0)*tor_s(0)
                                         + tor_r(1)*tor_s(1)
                                         + (tor(0)-Tor(0))*tor_rs(0)
                                         + (tor(1)-Tor(1))*tor_rs(1));
//              cout << "Kt" << endl << Kt << endl<< endl;

              //Total Stiffness
              K_Element(r,s) = K_Element(r,s) + (Km + Kb + Kt)*J1*Gauss_W(i);
						} // if (r <= s) {
					} // for (k1=0;k1<dof;k1++) {
				} // for (k=0;k<nos;k++) {
			} // for (j1=0;j1<dof;j1++)
		} // for (j=0;j<nos;j++)
	} // for (i=0;i<ngp;i++)
//    cout << "Residue Vector: "<< endl << F_Element << endl << endl;
//    cout << "Stiffness Matrix: "<< endl << K_Element << endl << endl;
  return make_tuple(K_Element,F_Element);
}
