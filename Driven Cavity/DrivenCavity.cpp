#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int N=20;
const int M=20;

//Deffinition of types
typedef double array1Dx[N];
typedef double array1Dy[M];
typedef double array2D[N][M];
typedef double staggeredx[N+1][M];
typedef double staggeredy[N][M+1];

//Deffinition of functions
void Geometry(int N, int M, double dx, double dy, array1Dx& x, array1Dy& y, array1Dy& Sx, array1Dx& Sy, array2D& V);
void ExportPoint(int N, int M, array1Dx x, array1Dy y, double xpoint, double ypoint, int& ipoint, int& jpoint);
void Property(int N, int M, double A, array2D& B);
void InitValueStagg(int N, int M, double A, staggeredx& Bx, staggeredy& By, double ALeft, double ARight, double ABot, double ATop);
void InitValue(int N, int M, double A, array2D& B);
void MassFlow(int N, int M, array1Dy Sx, array1Dx Sy, array2D rho, staggeredx vx, staggeredy vy, array2D& mw, array2D& me, array2D& ms, array2D& mn);
void ComputeRvx(int N, int M, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D V, staggeredx vx, double vxBot, double vxTop, double vxLeft, double vxRight, array2D mu, array2D mw, array2D me, array2D ms, array2D mn, staggeredx& Rvx);
void ComputeRvy(int N, int M, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D V, staggeredy vy, double vyBot, double vyTop, double vyLeft, double vyRight, array2D mu, array2D mw, array2D me, array2D ms, array2D mn, staggeredy& Rvy);
void ComputevPx(int N, int M, double dt, array2D rho, staggeredx vx, staggeredx Rvx, staggeredx Rvxm1, staggeredx& vPx);
void ComputevPy(int N, int M, double dt, array2D rho, staggeredy vy, staggeredy Rvy, staggeredy Rvym1, staggeredy& vPy);
void Divergence(int N, int M, double dx, double dy, staggeredx Ax, staggeredy Ay, array2D& DeltaA);
void CoefficientaWESN(int N, int M, array1Dy Sx, array1Dx Sy, double dx, double dy, double k, array2D& aW, array2D& aE, array2D& aS, array2D& aN);
void BoundaryAdiab(int N, int M, array2D& aW, array2D& aE, array2D& aS, array2D& aN);
void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D& aP);
void CoefficientbP(int N, int M, double dt, array2D V, array2D S, array2D rho, array2D& bP);
void Gradientx(int N, int M, double dx, array2D A, staggeredx& gradxA);
void Gradienty(int N, int M, double dy, array2D A, staggeredy& gradyA);
void Computevp1(int N, int M, double dt, array2D rho, staggeredx vPx, staggeredy vPy, staggeredx gradxP, staggeredy gradyP, staggeredx& vxp1, staggeredy& vyp1);
void NextValues(int N, int M, array2D Anew, array2D& Aold);
void NextValuesStagg(int N, int M, staggeredx Axnew, staggeredy Aynew, staggeredx& Axold, staggeredy& Ayold);
void ExportDataVX(int M, int ipoint, array1Dy y, staggeredx vxp1);
void ExportDataVY(int N, int jpoint, array1Dx x, staggeredy vyp1);
double TimeStep(double dx, double dy, double density, double difusivity, double vRef);
double MidPoint(double A, double B);
double MidPoint2D(double A1, double A2, double B1, double B2);
double Equation(double PhiW, double aW, double PhiE, double aE, double PhiS, double aS, double PhiN, double aN, double bP, double aP);

int main(){
	cout<<"PROGRAM INITIATIED"<<endl<<endl;
	
	//Geometry
	cout<<"Calculing geometry...";
	const double L=1;
	const double H=1;
	const double dx=L/double(N);
	const double dy=H/double(M);
	
	array1Dx x;
	array1Dy y;
	array1Dy Sx;
	array1Dx Sy;
	array2D V;
	
	Geometry(N, M, dx, dy, x, y, Sx, Sy, V);
	cout<<"Done"<<endl;
	
	//Export points
	double xpoint=0.5;
	double ypoint=0.5;
	int ipoint;
	int jpoint;
	
	ExportPoint(N, M, x, y, xpoint, ypoint, ipoint, jpoint);
	
	//Properties
	cout<<"Calculing properties...";
	double Re=1000;
	double vRef=1;
	double density=1;
	double difusivity=density*vRef*L/Re;
	
	array2D rho;
	array2D mu;
	
	Property(N, M, density, rho);
	Property(N, M, difusivity, mu);
	cout<<"Done"<<endl;
	
	//Time step
	double TimeEnd=50000;
	double dt;
	dt=TimeStep(dx, dy, density, difusivity, vRef);
	
	//Boundary values
	double vxLeft=0;
	double vxRight=0;
	double vxBot=0;
	double vxTop=1;
	double vyLeft=0;
	double vyRight=0;
	double vyBot=0;
	double vyTop=0;
	
	//Initial values
	cout<<"Initializating variables...";
	double vInit=0;
	
	staggeredx vx;
	staggeredy vy;
	staggeredx vxp1;
	staggeredy vyp1;
	
	InitValueStagg(N, M, vInit, vx, vy, vxLeft, vxRight, vyBot, vyTop);
	InitValueStagg(N, M, vInit, vxp1, vyp1, vxLeft, vxRight, vyBot, vyTop);
	
	double RvInit=0;
	
	staggeredx Rvxm1;
	staggeredy Rvym1;
	
	InitValueStagg(N, M, RvInit, Rvxm1, Rvym1, RvInit, RvInit, RvInit, RvInit);
	
	double PInit=0;
	
	array2D P;
	array2D PEst;
	
	InitValue(N, M, PInit, P);
	InitValue(N, M, PInit, PEst);
	cout<<"Done"<<endl;
	
	//Variables
	array2D mw;
	array2D me;
	array2D ms;
	array2D mn;	
	staggeredx Rvx;
	staggeredy Rvy;
	
	staggeredx vPx;
	staggeredy vPy;
	
	array2D divvP;
	
	array2D aP;
	array2D aW;
	array2D aE;
	array2D aS;
	array2D aN;
	array2D bP;
	double PW;
	double PE;
	double PS;
	double PN;
	staggeredx gradxP;
	staggeredy gradyP;
	
	//Additional data
	double delta=1e-6;
	double maxim=1;
	double Time=0;
	
	//Velocity field calculation
	cout<<"Calculatin velocity field over time...";
	while(Time<TimeEnd){
		cout<<Time<<endl;
		maxim=1;
		Time=Time+dt;
		MassFlow(N, M, Sx, Sy, rho, vx, vy, mw, me, ms, mn);
		ComputeRvx(N, M, dx, dy, Sx, Sy, V, vx, vxBot, vxTop, vxLeft, vxRight, mu, mw, me, ms, mn, Rvx);
		ComputeRvy(N, M, dx, dy, Sx, Sy, V, vy, vyBot, vyTop, vyLeft, vyRight, mu, mw, me, ms, mn, Rvy);
		ComputevPx(N, M, dt, rho, vx, Rvx, Rvxm1, vPx);
		ComputevPy(N, M, dt, rho, vy, Rvy, Rvym1, vPy);
		Divergence(N, M, dx, dy, vPx, vPy, divvP);
		CoefficientaWESN(N, M, Sx, Sy, dx, dy, 1, aW, aE, aS, aN);
		CoefficientaP(N, M, aW, aE, aS, aN, aP);
		CoefficientbP(N, M, dt, V, divvP, rho, bP);
		BoundaryAdiab(N, M, aW, aE, aS, aN);
		while(maxim>delta){
			maxim=0;
			for(int i=0; i<N; i++){
				for(int j=0; j<M; j++){
					PW=P[i-1][j];
					PE=P[i+1][j];
					PS=P[i][j-1];
					PN=P[i][j+1];
					if(i==0){
						PW=0;
					}
					if(i==N-1){
						PE=0;
					}
					if(j==0){
						PS=0;
					}
					if(j==M-1){
						PN=0;
					}
					P[i][j]=Equation(PW, aW[i][j], PE, aE[i][j], PS, aS[i][j], PN, aN[i][j], bP[i][j], aP[i][j]);
					if(fabs(P[i][j]-PEst[i][j])>maxim){
						maxim=fabs(P[i][j]-PEst[i][j]);
					}
					PEst[i][j]=P[i][j];
				}
			}
		}
		Gradientx(N, M, dx, P, gradxP);
		Gradienty(N, M, dy, P, gradyP);
		Computevp1(N, M, dt, rho, vPx, vPy, gradxP, gradyP, vxp1, vyp1);
		NextValuesStagg(N, M, Rvx, Rvy, Rvxm1, Rvym1);
		NextValuesStagg(N, M, vxp1, vyp1, vx, vy);
	}
	cout<<"Done"<<endl;
	
	//Data export
	cout<<"Exporting data...";
	ExportDataVX(M, ipoint, y, vxp1);
	ExportDataVY(N, jpoint, x, vyp1);
	cout<<"Done"<<endl;
	
	cout<<"PROGRAM FINISHED"<<endl;
	
	return 0;
}

void Geometry(int N, int M, double dx, double dy, array1Dx& x, array1Dy& y, array1Dy& Sx, array1Dx& Sy, array2D& V){
	for(int i=0; i<N; i++){
		x[i]=dx/2+i*dx;
		Sy[i]=dx;
	}
	for(int j=0; j<M; j++){
		y[j]=dy/2+j*dy;
		Sx[j]=dy;
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			V[i][j]=dx*dy;
		}
	}
}

void ExportPoint(int N, int M, array1Dx x, array1Dy y, double xpoint, double ypoint, int& ipoint, int& jpoint){
	for(int i=0; i<N; i++){
		if(x[i]<=xpoint&&x[i+1]>xpoint){
			ipoint=i;
		}
	}
	for(int j=0; j<M; j++){
		if(y[j]<=ypoint&&y[j+1]>ypoint){
			jpoint=j;
		}
	}
}

void Property(int N, int M, double A, array2D& B){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			B[i][j]=A;
		}
	}
}

void InitValueStagg(int N, int M, double A, staggeredx& Bx, staggeredy& By, double ALeft, double ARight, double ABot, double ATop){
	for(int i=0; i<N+1; i++){
		for(int j=0; j<M; j++){
			Bx[i][j]=A;
			if(i==0){
				Bx[i][j]=ALeft;
			}
			if(i==N){
				Bx[i][j]=ARight;
			}
		}
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<M+1; j++){
			By[i][j]=A;
			if(j==0){
				By[i][j]=ABot;
			}
			if(j==M){
				By[i][j]=ATop;
			}
		}
	}
}

void InitValue(int N, int M, double A, array2D& B){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			B[i][j]=A;
		}
	}
}

void MassFlow(int N, int M, array1Dy Sx, array1Dx Sy, array2D rho, staggeredx vx, staggeredy vy, array2D& mw, array2D& me, array2D& ms, array2D& mn){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			mw[i][j]=rho[i][j]*Sx[j]*vx[i][j];
			me[i][j]=rho[i][j]*Sx[j]*vx[i+1][j];
			ms[i][j]=rho[i][j]*Sy[i]*vy[i][j];
			mn[i][j]=rho[i][j]*Sy[i]*vy[i][j+1];
		}
	}
}

void ComputeRvx(int N, int M, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D V, staggeredx vx, double vxBot, double vxTop, double vxLeft, double vxRight, array2D mu, array2D mw, array2D me, array2D ms, array2D mn, staggeredx& Rvx){
	double vwx;
	double vex;
	double vsx;
	double vnx;
	double vPx;
	double vWx;
	double vEx;
	double vSx;
	double vNx;
	double Dw;
	double De;
	double Ds;
	double Dn;
	double Mfloww;
	double Mflowe;
	double Mflows;
	double Mflown;
	double mup=mu[1][1];
	for(int i=0; i<N+1; i++){
		for(int j=0; j<M; j++){
			vwx=MidPoint(vx[i][j], vx[i-1][j]);
			vex=MidPoint(vx[i][j], vx[i+1][j]);
			vsx=MidPoint(vx[i][j], vx[i][j-1]);
			vnx=MidPoint(vx[i][j], vx[i][j+1]);
			vPx=vx[i][j];
			vWx=vx[i-1][j];
			vEx=vx[i+1][j];
			vSx=vx[i][j-1];
			vNx=vx[i][j+1];
			Dw=mup*Sx[j]/dx;
			De=mup*Sx[j]/dx;
			Ds=mup*Sy[i]/dy;
			Dn=mup*Sy[i]/dy;
			Mfloww=MidPoint(mw[i-1][j], me[i-1][j]);
			Mflowe=MidPoint(mw[i][j], me[i][j]);
			Mflows=MidPoint(ms[i-1][j], ms[i][j]);
			Mflown=MidPoint(mn[i-1][j], mn[i][j]);
			if(i==0){
				vPx=vxLeft;
				vwx=0;
				vWx=0;
				Dw=0;
				Ds=mup*(Sy[i]/2)/dy;
				Dn=mup*(Sy[i]/2)/dy;
				Mfloww=0;
				Mflows=ms[i][j]/2;
				Mflown=mn[i][j]/2;
			}
			if(i==N){
				vPx=vxRight;
				vex=0;
				vEx=0;
				De=0;
				Ds=mup*(Sy[i-1]/2)/dy;
				Dn=mup*(Sy[i-1]/2)/dy;
				Mflowe=0;
				Mflows=ms[i][j]/2;
				Mflown=mn[i][j]/2;
			}
			if(j==0){
				vsx=vxBot;
				vSx=vxBot;
				Ds=mup*Sy[i]/(dy/2);
				if(i==0){
					Ds=mup*(Sy[i]/2)/(dy/2);
				}
				if(i==N){
					Ds=mup*(Sy[i-1]/2)/(dy/2);
				}
			}
			if(j==M-1){
				vnx=vxTop;
				vNx=vxTop;
				Dn=mup*Sy[i]/(dy/2);
				if(i==0){
					Dn=mup*(Sy[i]/2)/(dy/2);
				}
				if(i==N){
					Dn=mup*(Sy[i-1]/2)/(dy/2);
				}
			}
			Rvx[i][j]=Mfloww*vwx-Mflowe*vex+Mflows*vsx-Mflown*vnx+Dw*(vWx-vPx)+De*(vEx-vPx)+Ds*(vSx-vPx)+Dn*(vNx-vPx);
//			Rvx[i][j]=Rvx[i][j]/V[i][j];
//			if(i==0){
//				Rvx[i][j]=Mfloww*vwx-Mflowe*vex+Mflows*vsx-Mflown*vnx+Dw*(vWx-vPx)+De*(vEx-vPx)+Ds*(vSx-vPx)+Dn*(vNx-vPx);
//				Rvx[i][j]=Rvx[i][j]/(V[i][j]/2);
//			}
//			if(i==N){
//				Rvx[i][j]=Mfloww*vwx-Mflowe*vex+Mflows*vsx-Mflown*vnx+Dw*(vWx-vPx)+De*(vEx-vPx)+Ds*(vSx-vPx)+Dn*(vNx-vPx);
//				Rvx[i][j]=Rvx[i][j]/(V[i-1][j]/2);
//			}
		}
	}
}

void ComputeRvy(int N, int M, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D V, staggeredy vy, double vyBot, double vyTop, double vyLeft, double vyRight, array2D mu, array2D mw, array2D me, array2D ms, array2D mn, staggeredy& Rvy){
	double vwy;
	double vey;
	double vsy;
	double vny;
	double vPy;
	double vWy;
	double vEy;
	double vSy;
	double vNy;
	double Dw;
	double De;
	double Ds;
	double Dn;
	double Mfloww;
	double Mflowe;
	double Mflows;
	double Mflown;
	double mup=mu[1][1];
	for(int i=0; i<N; i++){
		for(int j=0; j<M+1; j++){
			vwy=MidPoint(vy[i][j], vy[i-1][j]);
			vey=MidPoint(vy[i][j], vy[i+1][j]);
			vsy=MidPoint(vy[i][j], vy[i][j-1]);
			vny=MidPoint(vy[i][j], vy[i][j+1]);
			vPy=vy[i][j];
			vWy=vy[i-1][j];
			vEy=vy[i+1][j];
			vSy=vy[i][j-1];
			vNy=vy[i][j+1];
			Dw=mup*Sx[j]/dx;
			De=mup*Sx[j]/dx;
			Ds=mup*Sy[i]/dy;
			Dn=mup*Sy[i]/dy;
			Mfloww=MidPoint(mw[i][j], mw[i][j-1]);
			Mflowe=MidPoint(me[i][j], me[i][j-1]);
			Mflows=MidPoint(ms[i][j-1], mn[i][j-1]);
			Mflown=MidPoint(ms[i][j], mn[i][j]);
			if(i==0){
				vwy=vyLeft;
				vWy=vyLeft;
				Dw=mup*Sx[j]/(dx/2);
			}
			if(i==N-1){
				vey=vyRight;
				vEy=vyRight;
				De=mup*Sx[j]/(dx/2);
			}
			if(j==0){
				vPy=vyBot;
				vsy=0;
				vSy=0;
				Dw=mup*(Sx[j]/2)/dx;
				De=mup*(Sx[j]/2)/dx;
				Ds=0;
				Mflows=0;
				Mfloww=mw[i][j]/2;
				Mflowe=me[i][j]/2;
				if(i==0){
					Dw=mup*(Sx[j]/2)/(dx/2);
				}
				if(i==N-1){
					De=mup*(Sx[j]/2)/(dx/2);
				}
			}
			if(j==M){
				vPy=vyRight;
				vny=0;
				vNy=0;
				Dw=mup*(Sx[j-1]/2)/dx;
				De=mup*(Sx[j-1]/2)/dx;
				Dn=0;
				Mflown=0;
				Mfloww=mw[i][j]/2;
				Mflowe=me[i][j]/2;
				if(i==0){
					Dw=mup*(Sx[j-1]/2)/(dx/2);
				}
				if(i==N-1){
					De=mup*(Sx[j-1]/2)/(dx/2);
				}
				
			}
			Rvy[i][j]=Mfloww*vwy-Mflowe*vey+Mflows*vsy-Mflown*vny+Dw*(vWy-vPy)+De*(vEy-vPy)+Ds*(vSy-vPy)+Dn*(vNy-vPy);
//			Rvy[i][j]=Rvy[i][j]/V[i][j];
//			if(j==0){
//				Rvy[i][j]=Mfloww*vwy-Mflowe*vey+Mflows*vsy-Mflown*vny+Dw*(vWy-vPy)+De*(vEy-vPy)+Ds*(vSy-vPy)+Dn*(vNy-vPy);
//				Rvy[i][j]=Rvy[i][j]/(V[i][j]/2);
//			}
//			if(j==M){
//				Rvy[i][j]=Mfloww*vwy-Mflowe*vey+Mflows*vsy-Mflown*vny+Dw*(vWy-vPy)+De*(vEy-vPy)+Ds*(vSy-vPy)+Dn*(vNy-vPy);
//				Rvy[i][j]=Rvy[i][j]/(V[i][j-1]/2);
//			}
		}
	}
}

void ComputevPx(int N, int M, double dt, array2D rho, staggeredx vx, staggeredx Rvx, staggeredx Rvxm1, staggeredx& vPx){
	double Rvxn;
	double Rvxnm1;
	double rhon;
	for(int i=0; i<N+1; i++){
		for(int j=0; j<M; j++){
			rhon=MidPoint(rho[i-1][j], rho[i][j]);
			Rvxn=Rvx[i][j];
			Rvxnm1=Rvxm1[i][j];
			if(i==0){
				rhon=rho[i][j];
			}
			if(i==N){
				rhon=rho[i-1][j];
			}
			vPx[i][j]=vx[i][j]+(dt/rhon)*((3*Rvxn/2)-(Rvxnm1/2));
		}
	}
}

void ComputevPy(int N, int M, double dt, array2D rho, staggeredy vy, staggeredy Rvy, staggeredy Rvym1, staggeredy& vPy){
	double Rvyn;
	double Rvynm1;
	double rhon;
	for(int i=0; i<N; i++){
		for(int j=0; j<M+1; j++){
			rhon=MidPoint(rho[i][j-1], rho[i][j]);
			Rvyn=Rvy[i][j];
			Rvynm1=Rvym1[i][j];
			if(j==0){
				rhon=rho[i][j];
			}
			if(j==M){
				rhon=rho[i][j-1];
			}
			vPy[i][j]=vy[i][j]+(dt/rhon)*((3*Rvyn/2)-(Rvynm1/2));
		}
	}
}

void Divergence(int N, int M, double dx, double dy, staggeredx Ax, staggeredy Ay, array2D& DeltaA){
	for(int i=0; i<N; i++){
		for(int j=0; j<M;j++){
			DeltaA[i][j]=(Ax[i+1][j]-Ax[i][j])/dx+(Ay[i][j+1]-Ay[i][j])/dy;
		}
	}
}

void CoefficientaWESN(int N, int M, array1Dy Sx, array1Dx Sy, double dx, double dy, double k, array2D& aW, array2D& aE, array2D& aS, array2D& aN){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			aW[i][j]=k*Sx[j]/dx;
			aE[i][j]=k*Sx[j]/dx;
			aS[i][j]=k*Sy[i]/dy;
			aN[i][j]=k*Sy[i]/dy;
			if(i==0){
				aW[i][j]=k*Sx[j]/(dx/2);
			}
			if(i==N-1){
				aE[i][j]=k*Sx[j]/(dx/2);
			}
			if(j==0){
				aS[i][j]=k*Sy[i]/(dy/2);
			}
			if(j==M-1){
				aN[i][j]=k*Sy[i]/(dy/2);
			}
		}
	}
}

void BoundaryAdiab(int N, int M, array2D& aW, array2D& aE, array2D& aS, array2D& aN){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			if(i==0){
				aW[i][j]=0;
			}
			if(i==N-1){
				aE[i][j]=0;
			}
			if(j==0){
				aS[i][j]=0;
			}
			if(j==M-1){
				aN[i][j]=0;
			}
		}
	}
}

void Gradientx(int N, int M, double dx, array2D A, staggeredx& gradxA){
	for(int i=0; i<N+1; i++){
		for(int j=0; j<M; j++){
			gradxA[i][j]=(A[i][j]-A[i-1][j])/dx;
			if(i==0){
				gradxA[i][j]=0;
			}
			if(i==N){
				gradxA[i][j]=0;
			}
		}
	}
}

void Gradienty(int N, int M, double dy, array2D A, staggeredy& gradyA){
	for(int i=0; i<N; i++){
		for(int j=0; j<M+1; j++){
			gradyA[i][j]=(A[i][j]-A[i][j-1])/dy;
			if(j==0){
				gradyA[i][j]=0;
			}
			if(j==M){
				gradyA[i][j]=0;
			}
		}
	}
}

void Computevp1(int N, int M, double dt, array2D rho, staggeredx vPx, staggeredy vPy, staggeredx gradxP, staggeredy gradyP, staggeredx& vxp1, staggeredy& vyp1){
	for(int i=1; i<N; i++){
		for(int j=0; j<M; j++){
			vxp1[i][j]=vPx[i][j]-dt*gradxP[i][j]/((rho[i][j]+rho[i-1][j])/2);
		}
	}
	for(int i=0; i<N; i++){
		for(int j=1; j<M; j++){
			vyp1[i][j]=vPy[i][j]-dt*gradyP[i][j]/((rho[i][j]+rho[i-1][j])/2);
		}
	}
}

void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D& aP){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			aP[i][j]=aW[i][j]+aE[i][j]+aS[i][j]+aN[i][j];
		}
	}
}

void CoefficientbP(int N, int M, double dt, array2D V, array2D S, array2D rho, array2D& bP){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			bP[i][j]=-rho[i][j]*S[i][j]*V[i][j]/dt;
		}
	}
}

void NextValues(int N, int M, array2D Anew, array2D& Aold){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			Aold[i][j]=Anew[i][j];
		}
	}
}

void NextValuesStagg(int N, int M, staggeredx Axnew, staggeredy Aynew, staggeredx& Axold, staggeredy& Ayold){
	for(int i=0; i<N+1; i++){
		for(int j=0; j<M; j++){
			Axold[i][j]=Axnew[i][j];
		}
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<M+1; j++){
			Ayold[i][j]=Aynew[i][j];
		}
	}
}

void ExportDataVX(int M, int ipoint, array1Dy y, staggeredx vxp1){
	ofstream file("VX.dat");
	for(int j=0; j<M; j++){
		file<<y[j]<<"	"<<vxp1[ipoint][j]<<endl;
	}
	file.close();
}

void ExportDataVY(int N, int jpoint, array1Dx x, staggeredy vyp1){
	ofstream file("VY.dat");
	for(int i=0; i<N; i++){
		file<<x[i]<<"	"<<vyp1[i][jpoint]<<endl;
	}
	file.close();
}

double TimeStep(double dx, double dy, double density, double difusivity, double vRef){
	double d;
	double dtc;
	double dtd;
	if(dy<dx){
		d=dy/2;
	}
	else{
		d=dx/2;
	}
	dtc=0.35*d/vRef;
	dtd=0.20*density*d*d/difusivity;
	if(dtc<dtd){
		return dtc;
	}
	else{
		return dtd;
	}
}

double MidPoint(double A, double B){
	double result=(A+B)/2;
	return result;
}

double MidPoint2D(double A1, double A2, double B1, double B2){
	double result1=(A1+A2)/2;
	double result2=(B1+B2)/2;
	double result=(result1+result2)/2;
	return result;
}

double Equation(double PhiW, double aW, double PhiE, double aE, double PhiS, double aS, double PhiN, double aN, double bP, double aP){
	double result=(aW*PhiW+aE*PhiE+aS*PhiS+aN*PhiN+bP)/aP;
	return result;
}

