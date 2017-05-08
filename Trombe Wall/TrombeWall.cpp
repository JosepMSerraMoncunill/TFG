#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int N1=2;
const int N2=10;
const int N3=4;
const int Ng=5;
const int T=10;

typedef double array1D[N1+N2+N3];
typedef double array1Dg[Ng];

void GeometryS(int N1, int N2, int N3, double dx1, double dx2, double dx3, double H, double W, array1D& x, array1D& Sy, array1D& V);
void GeometryG(int N, double dx, double H, double W, array1Dg& x, array1Dg& Sy, array1Dg& V);
void PropertyS(int N1, int N2, int N3, double phi1, double phi2, double phi3, array1D& phi);
void PropertyG(int N, double phi, array1Dg& phig);
void InitvalueS(int N, double B, array1D& A);
void InitvalueG(int N, double B, array1Dg& A);

int main(){
	//Constants definition
	const double pi=3.141592654;
	
	//Temporal data
	double TimeEnd=10;
	double dt=TimeEnd/double(T);
	
	//Geometry of the solid
	double ew1=0.050;
	double ew2=0.240;
	double ew3=0.010;
	double ea=0.100;
	double eg=0.004;
	double H=1;
	double W=1;
	
	double dx1=ew1/double(N1);
	double dx2=ew2/double(N2);
	double dx3=ew3/double(N3);
	double Sx=H*W;
	
	array1D xs;
	array1D Sy;
	array1D V;
	
	GeometryS(N1, N2, N3, dx1, dx2, dx3, H, W, xs, Sy, V);
	
	//Geometry of the glass
	double dxg=eg/double(Ng);
	double Sxg=H*W;
	
	array1Dg xg;
	array1Dg Syg;
	array1Dg Vg;
	
	GeometryG(Ng, dxg, H, W, xg, Syg, Vg);
	
	//Properties of the solid
	double lambda1=1;
	double lambda2=1;
	double lambda3=1;
	double rho1=1;
	double rho2=1;
	double rho3=1;
	
	array1D lambda;
	array1D rho;
	
	PropertyS(N1, N1, N3, lambda1, lambda2, lambda3, lambda);
	PropertyS(N1, N1, N3, rho1, rho2, rho3, rho);
	
	//Properties of the glass
	double lambda4=1;
	double rho4=1;
	
	array1Dg lambdag;
	array1Dg rhog;
	
	PropertyG(Ng, lambda4, lambdag);
	PropertyG(Ng, rho4, rhog);
	
	//Initial values of the solid
	double TInit=0;
	
	array1D T0;
	array1D T1;
	array1D TEst;
	
	InitvalueS(N1+N2+N3, TInit, T0);
	InitvalueS(N1+N2+N3, TInit, T1);
	InitvalueS(N1+N2+N3, TInit, TEst);
	
	//Initial values of the glass
	array1Dg T0g;
	array1Dg T1g;
	array1Dg TEstg;
	
	InitvalueG(Ng, TInit, T0g);
	InitvalueG(Ng, TInit, T1g);
	InitvalueG(Ng, TInit, TEstg);
	
	//Solar Radiation
	double IsDI=100;
	double thetaDI=30;
	double TsDI=0.7;
	double AsDI=0.2;
	double RsDI=0.1;
	
	double qradDI1=-TsDI*IsDI*cos(thetaDI*pi/180);
	double qradDI2=TsDI*IsDI*cos(thetaDI*pi/180);
	double qradDI3=RsDI*IsDI*cos(thetaDI*pi/180)-IsDI*cos(thetaDI*pi/180);
	
	double IsDF=100;
	double TsDF=0.7;
	double AsDF=0.2;
	double RsDF=0.1;
	double rsDF1=0.1;
	
	double qradDF1;
	double qradDF2;
	double qradDF3;
	
	return 0;
}

void GeometryS(int N1, int N2, int N3, double dx1, double dx2, double dx3, double H, double W, array1D& x, array1D& Sy, array1D& V){
	for(int i=0; i<N1+N2+N3; i++){
		x[i]=dx1/2+i*dx1;
		Sy[i]=dx1*H;
		V[i]=dx1*H*W;
	}
	for(int i=0; i<N2; i++){
		x[N1-1+i]=x[N1-1]+dx1/2+dx2/2+i*dx2;
		Sy[N1-1+i]=dx2*H;
		V[N1-1+1]=dx2*H*W;
	}
	for(int i=0; i<N3; i++){
		x[N1+N2-1+i]=x[N1+N2-1]+dx2/2+dx3/3+i*dx3;
		Sy[N1+N2-1+i]=dx3*H;
		V[N1+N2-1+i]=dx3*H*W;
	}
}

void GeometryG(int N, double dx, double H, double W, array1Dg& x, array1Dg& Sy, array1Dg& V){
	for(int i=0; i<N; i++){
		x[i]=dx/2+i*dx;
		Sy[i]=dx*H;
		V[i]=dx*H*W;
	}
}

void PropertyS(int N1, int N2, int N3, double phi1, double phi2, double phi3, array1D& phi){
	for(int i=0; i<N1+N2+N3; i++){
		if(i<N1){
			phi[i]=phi1;
		}
		if(i>=N1&&i<N2){
			phi[i]=phi2;
		}
		if(i>=N2){
			phi[i]=phi3;
		}
	}
}

void PropertyG(int N, double phi, array1Dg& phig){
	for(int i=0; i<N; i++){
		phig[i]=phi;
	}
}

void InitvalueS(int N, double B, array1D& A){
	for(int i=0; i<N; i++){
		A[i]=B;
	}
}

void InitvalueG(int N, double B, array1Dg& A){
	for(int i=0; i<N; i++){
		A[i]=B;
	}
}
