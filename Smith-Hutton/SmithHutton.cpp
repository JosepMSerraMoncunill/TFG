#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int N=160;
const int M=100;
const int Np1=N+1;
const int Mp1=M+1;
const int Time=10000;
const double TimeEnd=10000;

//Deffinition of types
typedef double array1Dx[N];
typedef double array1Dy[M];
typedef double array2D[N][M];

//Deffinition of functions
void Geometry(int N, int M, double dx, double dy, array1Dx& x, array1Dy& y, array1Dy& Sx, array1Dx& Sy, array2D& V);
void ExportPoint(int N, array1Dx x, double xpoint, int& ipoint);
void Properties(int N, int M, double density, double difusivity, array2D rho, array2D Gamma);
void InitPhi(int N, int M, double PhiInit, array2D& phi);
void MassFlow(int N, int M, double dx, double dy, array1Dx x, array1Dy y, array1Dy Sx, array1Dx Sy, array2D rho, array2D& mw, array2D& me, array2D& ms, array2D& mn);
void CoefficientWESN(int N, int M, array2D Gamma, double dx, double dy, array1Dx Sy, array1Dy Sx, array2D mw, array2D me, array2D ms, array2D mn, int scheme, array2D& aW, array2D& aE, array2D& aS, array2D& aN);
double AP(int scheme, double P);
void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D rho, array2D V, double dt, double Sp, array2D& aP);
void CoefficientbP(int i, int j, array2D rho, array2D V, array2D Phi0, double dt, double Sc, array2D& bP);
void PhiBoundary(double& Phi);
void PhiInlet(array1Dx x, int i, double& Phi);
void Equation(double PhiW, double PhiE, double PhiS, double PhiN, double aW, double aE, double aS, double aN, double aP, double bP, double& PhiP);
void ExportData(int M, array1Dx x, int ipoint, array2D Phi);
void ExportPlot(int N, int M, array1Dx x, array1Dy y, array2D Phi);
double max(double A, double B);

int main(){
	cout<<"PROGRAM INITIATIED"<<endl<<endl;
	
	//Geometry
	cout<<"Calculing geometry...";
	const double L=2;
	const double H=1;
	const double dx=L/double(N);
	const double dy=H/double(M);
	double dt=TimeEnd/double(Time);
	
	array1Dx x;
	array1Dy y;
	array1Dy Sx;
	array1Dx Sy;
	array2D V;
	
	Geometry(N, M, dx, dy, x, y, Sx, Sy, V);
	cout<<"Done"<<endl;
	
	//Export points
	double xpoint=0;
	int ipoint;
	ExportPoint(N, x, xpoint, ipoint);
	
	//Properties
	cout<<"Calculing properties...";
	const double density=1000;
	const double difusivity=1;
	const double Sp=0;
	const double Sc=0;
	
	array2D rho;
	array2D Gamma;
	
	Properties(N, M, density, difusivity, rho, Gamma);
	cout<<"Done"<<endl;
	
	//Mass flow
	cout<<"Calculing mass flows...";
	array2D mw;
	array2D me;
	array2D ms;
	array2D mn;
	
	MassFlow(N, M, dx, dy, x, y, Sx, Sy, rho, mw, me, ms, mn);
	cout<<"Done"<<endl;
	
	//Initialization of the variable
	cout<<"Assigning initial values...";
	const double PhiInit=0;
	
	array2D Phi0;
	array2D Phi;
	array2D PhiEst;
	
	InitPhi(N, M, PhiInit, Phi0);
	InitPhi(N, M, PhiInit, Phi);
	InitPhi(N, M, PhiInit, PhiEst);
	cout<<"Done"<<endl;
	
	//Coefficients
	cout<<"Calculing coefficients...";
	//Select a scheme from below:
	//	0 : Central difference
	//	1 : Upwind
	//	2 : Hybrid
	//	3 : Power law
	//	4 : Exponential
	int scheme=3;
	array2D aW;
	array2D aE;
	array2D aS;
	array2D aN;
	array2D aP;
	array2D bP;
	
	CoefficientWESN(N, M, Gamma, dx, dy, Sy, Sx, mw, me, ms, mn, scheme, aW, aE, aS, aN);
	CoefficientaP(N, M, aW, aE, aS, aN, rho, V, dt, Sp, aP);
	cout<<"Done"<<endl;
	
	double time=0;
	double delta=1e-5;
	double fr=1;
	double maximum;
	
	cout<<"Calculating Phi over time...";
	for(int t=0; t<Time; t++){
		maximum=1;
		while(maximum>delta){
			maximum=0;
			for(int i=0; i<N; i++){
				for(int j=0; j<M; j++){
					CoefficientbP(i, j, rho, V, Phi0, dt, Sc, bP);
					double PhiW=Phi[i-1][j];
					double PhiE=Phi[i+1][j];
					double PhiS=Phi[i][j-1];
					double PhiN=Phi[i][j+1];
					if(i==0){
						PhiBoundary(PhiW);
					}
					if(i==N-1){
						PhiBoundary(PhiE);
					}
					if(j==M-1){
						PhiBoundary(PhiN);
					}
					if(j==0&&i<=(N-1)/2){
						PhiInlet(x, i, PhiS);
					}
					if(j==0&&i>(N-1)/2){
						PhiS=Phi[i][j];
					}
					Equation(PhiW, PhiE, PhiS, PhiN, aW[i][j], aE[i][j], aS[i][j], aN[i][j], aP[i][j], bP[i][j], Phi[i][j]);
					
					if(fabs(Phi[i][j]-PhiEst[i][j])>maximum){
						maximum=fabs(Phi[i][j]-PhiEst[i][j]);
					}
					PhiEst[i][j]=PhiEst[i][j]+fr*(Phi[i][j]-PhiEst[i][j]);
				}
			}
		} //End of loop
		for(int i=0;i<N;i++){
			for(int j=0;j<M;j++){
				Phi0[i][j]=Phi[i][j];
			}
		}
		time=time+dt;
	}
	cout<<"Done"<<endl;
	
	cout<<"Exporting data...";
	ExportData(M, x, ipoint, Phi);
	cout<<"Done"<<endl;
	
	cout<<"Preparing data for plot...";
	ExportPlot(N, M, x, y, Phi);
	cout<<"Done"<<endl;
	
	cout<<endl<<"PROGRAM FINISHED"<<endl;
	
	return 0;
}

void Geometry(int N, int M, double dx, double dy, array1Dx& x, array1Dy& y, array1Dy& Sx, array1Dx& Sy, array2D& V){
	//Coordinates in x direction
	for(int i=0; i<N; i++){
		x[i]=-1+dx/2+i*dx;
	}
	//Coordinates in y direction
	for(int j=0; j<M; j++){
		y[j]=dy/2+j*dy;
	}
	//Surfaces normals to x direction
	for(int j=0; j<M; j++){
		Sx[j]=dy;
	}
	//Surfaces normals to y direction
	for(int i=0; i<N; i++){
		Sy[i]=dx;
	}
	//Volumes
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			V[i][j]=dx*dy;
		}
	}
}

void ExportPoint(int N, array1Dx x, double xpoint, int& ipoint){
	for(int i=0; i<N; i++){
		if(x[i]<=xpoint&&x[i+1]>xpoint){
			ipoint=i+1;
		}
	}
}

void Properties(int N, int M, double density, double difusivity, array2D rho, array2D Gamma){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			rho[i][j]=density;
			Gamma[i][j]=difusivity;
		}
	}
}

void MassFlow(int N, int M, double dx, double dy, array1Dx x, array1Dy y, array1Dy Sx, array1Dx Sy, array2D rho, array2D& mw, array2D& me, array2D& ms, array2D& mn){
	double xp1[N+1];
	for(int i=0; i<N+1; i++){
		xp1[i]=-1+dx*i;
	}
	double yp1[M+1];
	for(int j=0; j<M+1; j++){
		yp1[j]=dy*j;
	}
	double vw;
	double ve;
	double vs;
	double vn;
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			vw=2*y[j]*(1-xp1[i]*xp1[i]);
			ve=2*y[j]*(1-xp1[i+1]*xp1[i+1]);
			vs=-2*x[i]*(1-yp1[j]*yp1[j]);
			vn=-2*x[i]*(1-yp1[j+1]*yp1[j+1]);
			mw[i][j]=rho[i][j]*vw*Sx[j];
			me[i][j]=rho[i][j]*ve*Sx[j];
			ms[i][j]=rho[i][j]*vs*Sy[i];
			mn[i][j]=rho[i][j]*vn*Sy[i];
		}
	}
}

void InitPhi(int N, int M, double PhiInit, array2D& Phi){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			Phi[i][j]=PhiInit;
		}
	}
}

void CoefficientWESN(int N, int M, array2D Gamma, double dx, double dy, array1Dx Sy, array1Dy Sx, array2D mw, array2D me, array2D ms, array2D mn, int scheme, array2D& aW, array2D& aE, array2D& aS, array2D& aN){	
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			double Gammaw=dx/(0.5*dx/Gamma[i][j]+0.5*dx/Gamma[i-1][j]);
			double Gammae=dx/(0.5*dx/Gamma[i][j]+0.5*dx/Gamma[i+1][j]);
			double Gammas=dy/(0.5*dy/Gamma[i][j]+0.5*dy/Gamma[i][j-1]);
			double Gamman=dy/(0.5*dy/Gamma[i][j]+0.5*dy/Gamma[i][j+1]);
			double Dw=Gammaw*Sx[j]/dx;
			double De=Gammae*Sx[j]/dx;
			double Ds=Gammas*Sy[i]/dy;
			double Dn=Gamman*Sy[i]/dy;
			if(i==0){
				Gammaw=Gamma[i][j];
				Dw=Gammaw*Sx[j]/(0.5*dx);
			}
			if(i==N-1){
				Gammae=Gamma[i][j];
				De=Gammae*Sx[j]/(0.5*dx);
			}
			if(j==0){
				Gammas=Gamma[i][j];
				Ds=Gammas*Sy[i]/(0.5*dy);
			}
			if(j==M-1){
				Gamman=Gamma[i][j];
				Dn=Gamman*Sy[i]/(0.5*dy);
			}
			double Pw=mw[i][j]/Dw;
			double Pe=me[i][j]/De;
			double Ps=ms[i][j]/Ds;
			double Pn=mn[i][j]/Dn;
			double Aw=AP(scheme, Pw);
			double Ae=AP(scheme, Pe);
			double As=AP(scheme, Ps);
			double An=AP(scheme, Pn);
			aW[i][j]=Dw*Aw+max(mw[i][j],0);
			aE[i][j]=De*Ae+max(-me[i][j],0);
			aS[i][j]=Ds*As+max(ms[i][j],0);
			aN[i][j]=Dn*An+max(-mn[i][j],0);
		}
	}
}

double AP(int scheme, double P){
	if(scheme==0){
		return 1-0.5*fabs(P);
	}
	else if(scheme==1){
		return 1;
	}
	else if(scheme==2){
		return max(0,1-0.5*fabs(P));
	}
	else if(scheme==3){
		return max(0,pow(1-0.1*fabs(P),5));
	}
	else{
		return fabs(P)/(exp(fabs(P))-1);
	}
}

void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D rho, array2D V, double dt, double Sp, array2D& aP){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			aP[i][j]=aW[i][j]+aE[i][j]+aS[i][j]+aN[i][j]+rho[i][j]*V[i][j]/dt-Sp*V[i][j];
		}
	}
}

void CoefficientbP(int i, int j, array2D rho, array2D V, array2D Phi0, double dt, double Sc, array2D& bP){
	bP[i][j]=Sc*V[i][j]+Phi0[i][j]*rho[i][j]*V[i][j]/dt;
}

void PhiBoundary(double& Phi){
	double alpha=10;
	Phi=1-tanh(alpha);
}

void PhiInlet(array1Dx x, int i, double& Phi){
	double alpha=10;
	Phi=1+tanh(alpha*(2*x[i]+1));
}

void Equation(double PhiW, double PhiE, double PhiS, double PhiN, double aW, double aE, double aS, double aN, double aP, double bP, double& PhiP){
	PhiP=(PhiW*aW+PhiE*aE+PhiS*aS+PhiN*aN+bP)/aP;
}

void ExportData(int M, array1Dx x, int ipoint, array2D Phi){
	ofstream file("SERRA.dat");
	for(int i=ipoint; i<N; i++){
		file<<x[i]<<"	"<<(Phi[i][0]+Phi[i+1][0])/2<<endl;
	}
	file.close();
}

void ExportPlot(int N, int M, array1Dx x, array1Dy y, array2D Phi){
	ofstream file("DATA.dat");
	for(int j=0; j<M; j++){
		for(int i=0; i<N; i++){
			file<<x[i]<<"	"<<y[j]<<"	"<<Phi[i][j]<<endl;
		}
	}
	file.close();
}

double max(double A, double B){
	if(A>=B){
		return A;
	}
	else{
		return B;
	}
}

