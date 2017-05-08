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
void CoefficientWESN(int N, int M, array2D Gamma, double dx, double dy, array1Dx Sy, array1Dy Sx, array2D mw, array2D me, array2D ms, array2D mn, array2D& aW, array2D& aE, array2D& aS, array2D& aN);
void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D rho, array2D V, double dt, double Sp, array2D& aP);
void schemeUDS(int i, int j, array1Dx x, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, double& PhiwUDS, double& PhieUDS, double& PhisUDS, double& PhinUDS);
void schemeCDS(int i, int j, array1Dx x, array2D Phi, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS);
void schemeEDS(int i, int j, array1Dx x, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D Phi, array2D Gamma, array2D mw, array2D me, array2D ms, array2D mn, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS);
void schemeSUDS(int i, int j, double dx, double dy, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, array1Dx x, array1Dy y, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS);
void schemeQUICK(int i, int j, double dx, double dy, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, array1Dx x, array1Dy y, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS);
void CoefficientbP1(int i, int j, array2D mw, array2D me, array2D ms, array2D mn, double PhiwUDS, double PhieUDS, double PhisUDS, double PhinUDS, double PhiwHRS, double PhieHRS, double PhisHRS, double PhinHRS, array2D& bP);
void CoefficientbP2(int i, int j, array2D rho, array2D V, array2D Phi0, double dt, double Sp, array2D& bP);
void PhiBoundary(double& Phi);
void PhiInlet(array1Dx x, int i, double& Phi);
void Equation(double PhiW, double PhiE, double PhiS, double PhiN, double aW, double aE, double aS, double aN, double aP, double bP, double& PhiP);
void ExportData(int M, array1Dx x, int ipoint, array2D Phi);
void ExportPlot(int N, int M, array1Dx x, array1Dy y, array2D Phi);

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
	array2D aW;
	array2D aE;
	array2D aS;
	array2D aN;
	array2D aP;
	array2D bP;
	
	CoefficientWESN(N, M, Gamma, dx, dy, Sy, Sx, mw, me, ms, mn, aW, aE, aS, aN);
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
					double PhiwUDS;
					double PhieUDS;
					double PhisUDS;
					double PhinUDS;
					double PhiwHRS;
					double PhieHRS;
					double PhisHRS;
					double PhinHRS;
					schemeUDS(i, j, x, Phi, mw, me, ms, mn, PhiwUDS, PhieUDS, PhisUDS, PhinUDS);
					//choose 1 from below
					schemeUDS(i, j, x, Phi, mw, me, ms, mn, PhiwHRS, PhieHRS, PhisHRS, PhinHRS);
					//schemeCDS(i, j, x, Phi, PhiwHRS, PhieHRS, PhisHRS, PhinHRS);
					//schemeEDS(i, j, x, dx, dy, Sx, Sy, Phi, Gamma, mw, me, ms, mn, PhiwHRS, PhieHRS, PhisHRS, PhinHRS);
					//schemeSUDS(i, j, dx, dy, Phi, mw, me, ms, mn, x, y, PhiwHRS, PhieHRS, PhisHRS, PhinHRS);
					//schemeQUICK(i, j, dx, dy, Phi, mw, me, ms, mn, x, y, PhiwHRS, PhieHRS, PhisHRS, PhinHRS);
					CoefficientbP1( i, j, mw, me, ms, mn, PhiwUDS, PhieUDS, PhisUDS, PhinUDS, PhiwHRS, PhieHRS, PhisHRS, PhinHRS, bP);
					CoefficientbP2( i, j, rho, V, Phi0, dt, Sp, bP);
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
					if(j==0){
						PhiBoundary(PhiS);
					}
					if(j==M-1){
						PhiBoundary(PhiN);
					}
					if(i>=0&&i<=0.5*(N-1)&&j==0){
						PhiInlet(x, i, PhiS);
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

void CoefficientWESN(int N, int M, array2D Gamma, double dx, double dy, array1Dx Sy, array1Dy Sx, array2D mw, array2D me, array2D ms, array2D mn, array2D& aW, array2D& aE, array2D& aS, array2D& aN){	
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
			if(i==0){ //Left nodes
				Gammaw=Gamma[i][j];
				Dw=Gammaw*Sx[j]/(0.5*dx);
			}
			if(i==N-1){ //Right nodes
				Gammae=Gamma[i][j];
				De=Gammae*Sx[j]/(0.5*dx);
			}
			if(j==0){ //Bottom nodes
				Gammas=Gamma[i][j];
				Ds=Gammas*Sy[i]/(0.5*dy);
			}
			if(j==M-1){ //Top nodes
				Gamman=Gamma[i][j];
				Dn=Gamman*Sy[i]/(0.5*dy);
			}
			aW[i][j]=Dw+0.5*(mw[i][j]+fabs(mw[i][j]));
			aE[i][j]=De-0.5*(me[i][j]-fabs(me[i][j]));
			aS[i][j]=Ds+0.5*(ms[i][j]+fabs(ms[i][j]));
			aN[i][j]=Dn-0.5*(mn[i][j]-fabs(mn[i][j]));
		}
	}
}

void CoefficientaP(int N, int M, array2D aW, array2D aE, array2D aS, array2D aN, array2D rho, array2D V, double dt, double Sp, array2D& aP){
	for(int i=0; i<N; i++){
		for(int j=0; j<M; j++){
			aP[i][j]=aW[i][j]+aE[i][j]+aS[i][j]+aN[i][j]+rho[i][j]*V[i][j]/dt-Sp*V[i][j];
		}
	}
}

void schemeUDS(int i, int j, array1Dx x, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, double& PhiwUDS, double& PhieUDS, double& PhisUDS, double& PhinUDS){
	if(mw[i][j]>0){
		PhiwUDS=Phi[i-1][j];
	}
	else{
		PhiwUDS=Phi[i][j];
	}
	if(me[i][j]>0){
		PhieUDS=Phi[i][j];
	}
	else{
		PhieUDS=Phi[i+1][j];
	}
	if(ms[i][j]>0){
		PhisUDS=Phi[i][j-1];
	}
	else{
		PhisUDS=Phi[i][j];
	}
	if(mn[i][j]>0){
		PhinUDS=Phi[i][j];
	}
	else{
		PhinUDS=Phi[i][j+1];
	}
	if(i==0){
		PhiBoundary(PhiwUDS);
	}
	if(i==N-1){
		PhiBoundary(PhieUDS);
	}
	if(j==M-1){
		PhiBoundary(PhinUDS);
	}
	if(j==0&&i<=(N-1)/2){
		PhiInlet(x, i, PhisUDS);
	}
	if(j==0&&i>(N-1)/2){
		PhisUDS=Phi[i][j];
	}
}

void schemeCDS(int i, int j, array1Dx x, array2D Phi, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS){
	PhiwHRS=(Phi[i][j]+Phi[i-1][j])/2;
	PhieHRS=(Phi[i][j]+Phi[i+1][j])/2;
	PhisHRS=(Phi[i][j]+Phi[i][j-1])/2;
	PhinHRS=(Phi[i][j]+Phi[i][j+1])/2;
	if(i==0){
		PhiBoundary(PhiwHRS);
	}
	if(i==N-1){
		PhiBoundary(PhieHRS);
	}
	if(j==M-1){
		PhiBoundary(PhinHRS);
	}
	if(j==0&&i<=(N-1)/2){
		PhiInlet(x, i, PhisHRS);
	}
	if(j==0&&i>(N-1)/2){
		PhisHRS=Phi[i][j];
	}
}

void schemeEDS(int i, int j, array1Dx x, double dx, double dy, array1Dy Sx, array1Dx Sy, array2D Phi, array2D Gamma, array2D mw, array2D me, array2D ms, array2D mn, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS){
	double Gammaw=dx/(0.5*dx/Gamma[i][j]+0.5*dx/Gamma[i-1][j]);
	double Gammae=dx/(0.5*dx/Gamma[i][j]+0.5*dx/Gamma[i+1][j]);
	double Gammas=dy/(0.5*dy/Gamma[i][j]+0.5*dy/Gamma[i][j-1]);
	double Gamman=dy/(0.5*dy/Gamma[i][j]+0.5*dy/Gamma[i][j+1]);
	double Dw=Gammaw*Sx[j]/dx;
	double De=Gammae*Sx[j]/dx;
	double Ds=Gammas*Sy[i]/dy;
	double Dn=Gamman*Sy[i]/dy;
	if(i==0){ //Left nodes
		Gammaw=Gamma[i][j];
		Dw=Gammaw*Sx[j]/(0.5*dx);
	}
	if(i==N-1){ //Right nodes
		Gammae=Gamma[i][j];
		De=Gammae*Sx[j]/(0.5*dx);
	}
	if(j==0){ //Bottom nodes
		Gammas=Gamma[i][j];
		Ds=Gammas*Sy[i]/(0.5*dy);
	}
	if(j==M-1){ //Top nodes
		Gamman=Gamma[i][j];
		Dn=Gamman*Sy[i]/(0.5*dy);
	}
	double Pw=mw[i][j]/Dw;
	double Pe=me[i][j]/De;
	double Ps=ms[i][j]/Ds;
	double Pn=mn[i][j]/Dn;
	PhiwHRS=Phi[i-1][j]+(Phi[i][j]-Phi[i-1][j])*(exp(Pw/2)-1)/(exp(Pw)-1);
	PhieHRS=Phi[i][j]+(Phi[i+1][j]-Phi[i][j])*(exp(Pe/2)-1)/(exp(Pe)-1);
	PhisHRS=Phi[i][j-1]+(Phi[i][j]-Phi[i][j-1])*(exp(Ps/2)-1)/(exp(Ps)-1);
	PhinHRS=Phi[i][j]+(Phi[i][j+1]-Phi[i][j])*(exp(Pn/2)-1)/(exp(Pn)-1);
	if(i==0){
		PhiBoundary(PhiwHRS);
	}
	if(i==N-1){
		PhiBoundary(PhieHRS);
	}
	if(j==M-1){
		PhiBoundary(PhinHRS);
	}
	if(j==0&&i<=(N-1)/2){
		PhiInlet(x, i, PhisHRS);
	}
	if(j==0&&i>(N-1)/2){
		PhisHRS=Phi[i][j];
	}
}

void schemeSUDS(int i, int j, double dx, double dy, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, array1Dx x, array1Dy y, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS){
	double PhiC;
	double PhiU;
	double PhiD;
	double xC;
	double xU;
	double xD;
	double xF;
	double xAC;
	double xAF;
	double PhiAC;
	double PhiAF;
	if(mw>=0){
		PhiC=Phi[i-1][j];
		PhiU=Phi[i-2][j];
		PhiD=Phi[i][j];
		xC=x[i-1];
		xU=x[i-2];
		xD=x[i];
		xF=x[i]-dx/2;
		if(i==1){
			PhiBoundary(PhiU);
			xU=-1;
		}
	}
	else{
		PhiC=Phi[i][j];
		PhiU=Phi[i+1][j];
		PhiD=Phi[i-1][j];
		xC=x[i];
		xU=x[i+1];
		xD=x[i-1];
		xF=x[i]-dx/2;
		if(i==N-1){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=PhiAC*xAF/xAC;
	PhiwHRS=PhiAF*(PhiD-PhiU)+PhiU;
	if(me>=0){
		PhiC=Phi[i][j];
		PhiU=Phi[i-1][j];
		PhiD=Phi[i+1][j];
		xC=x[i];
		xU=x[i-1];
		xD=x[i+1];
		xF=x[i]+dx/2;
		if(i==0){
			PhiBoundary(PhiU);
			xU=-1;
		}
	}
	else{
		PhiC=Phi[i+1][j];
		PhiU=Phi[i+2][j];
		PhiD=Phi[i][j];
		xC=x[i+1];
		xU=x[i+2];
		xD=x[i];
		xF=x[i]+dx/2;
		if(i==N-2){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhieHRS=PhiAC*xAF/xAC;
	if(ms>=0){
		PhiC=Phi[i][j-1];
		PhiU=Phi[i][j-2];
		PhiD=Phi[i][j];
		xC=y[j-1];
		xU=y[j-2];
		xD=y[j];
		xF=y[j]-dy/2;
		if(j==1){
			if(i<=(N-1)/2){
				PhiInlet(x, i, PhiU);
			}
			else{
				PhiU=PhiC;
			}
			xU=0;
		}
	}
	else{
		PhiC=Phi[i][j];
		PhiU=Phi[i][j+1];
		PhiD=Phi[i][j-1];
		xC=y[j];
		xU=y[j+1];
		xD=y[j-1];
		xF=y[j]-dy/2;
		if(j==M-1){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhisHRS=PhiAC*xAF/xAC;
	if(mn>=0){
		PhiC=Phi[i][j];
		PhiU=Phi[i][j-1];
		PhiD=Phi[i][j+1];
		xC=y[j];
		xU=y[j-1];
		xD=y[j+1];
		xF=y[j]+dy/2;
		if(j==0){
			if(i<=(N-1)/2){
				PhiInlet(x, i, PhiU);
			}
			else{
				PhiU=PhiC;
			}
			PhiBoundary(PhiU);
			xU=0;
		}
	}
	else{
		PhiC=Phi[i][j+1];
		PhiU=Phi[i][j+2];
		PhiD=Phi[i][j];
		xC=y[j+1];
		xU=y[j+2];
		xD=y[j];
		xF=y[j]+dy/2;
		if(j==M-2){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhinHRS=PhiAC*xAF/xAC;
	if(i==0){
		PhiBoundary(PhiwHRS);
	}
	if(i==N-1){
		PhiBoundary(PhieHRS);
	}
	if(j==M-1){
		PhiBoundary(PhinHRS);
	}
	if(i<=0.5*(N-1)&&j==0){
		PhiInlet(x, i, PhisHRS);
	}
	if(i>0.5*(N-1)&&j==0){
		PhisHRS=Phi[i][j];
	}
}

void schemeQUICK(int i, int j, double dx, double dy, array2D Phi, array2D mw, array2D me, array2D ms, array2D mn, array1Dx x, array1Dy y, double& PhiwHRS, double& PhieHRS, double& PhisHRS, double& PhinHRS){
	double PhiC;
	double PhiU;
	double PhiD;
	double xC;
	double xU;
	double xD;
	double xF;
	double xAC;
	double xAF;
	double PhiAC;
	double PhiAF;
	if(mw>=0){
		PhiC=Phi[i-1][j];
		PhiU=Phi[i-2][j];
		PhiD=Phi[i][j];
		xC=x[i-1];
		xU=x[i-2];
		xD=x[i];
		xF=x[i]-dx/2;
		if(i==1){
			PhiBoundary(PhiU);
			xU=-1;
		}
	}
	else{
		PhiC=Phi[i][j];
		PhiU=Phi[i+1][j];
		PhiD=Phi[i-1][j];
		xC=x[i];
		xU=x[i+1];
		xD=x[i-1];
		xF=x[i]-dx/2;
		if(i==N-1){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhiwHRS=PhiAF*(PhiD-PhiU)+PhiU;
	if(me>=0){
		PhiC=Phi[i][j];
		PhiU=Phi[i-1][j];
		PhiD=Phi[i+1][j];
		xC=x[i];
		xU=x[i-1];
		xD=x[i+1];
		xF=x[i]+dx/2;
		if(i==0){
			PhiBoundary(PhiU);
			xU=-1;
		}
	}
	else{
		PhiC=Phi[i+1][j];
		PhiU=Phi[i+2][j];
		PhiD=Phi[i][j];
		xC=x[i+1];
		xU=x[i+2];
		xD=x[i];
		xF=x[i]+dx/2;
		if(i==N-2){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhieHRS=PhiAF*(PhiD-PhiU)+PhiU;
	if(ms>=0){
		PhiC=Phi[i][j-1];
		PhiU=Phi[i][j-2];
		PhiD=Phi[i][j];
		xC=y[j-1];
		xU=y[j-2];
		xD=y[j];
		xF=y[j]-dy/2;
		if(j==1){
			if(i<=(N-1)/2){
				PhiInlet(x, i, PhiU);
			}
			else{
				PhiU=PhiC;
			}
			xU=0;
		}
	}
	else{
		PhiC=Phi[i][j];
		PhiU=Phi[i][j+1];
		PhiD=Phi[i][j-1];
		xC=y[j];
		xU=y[j+1];
		xD=y[j-1];
		xF=y[j]-dy/2;
		if(j==M-1){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhisHRS=PhiAF*(PhiD-PhiU)+PhiU;
	if(mn>=0){
		PhiC=Phi[i][j];
		PhiU=Phi[i][j-1];
		PhiD=Phi[i][j+1];
		xC=y[j];
		xU=y[j-1];
		xD=y[j+1];
		xF=y[j]+dy/2;
		if(j==0){
			if(i<=(N-1)/2){
				PhiInlet(x, i, PhiU);
			}
			else{
				PhiU=PhiC;
			}
			PhiBoundary(PhiU);
			xU=0;
		}
	}
	else{
		PhiC=Phi[i][j+1];
		PhiU=Phi[i][j+2];
		PhiD=Phi[i][j];
		xC=y[j+1];
		xU=y[j+2];
		xD=y[j];
		xF=y[j]+dy/2;
		if(j==M-2){
			PhiBoundary(PhiU);
			xU=1;
		}
	}
	xAC=(xC-xU)/(xD-xU);
	xAF=(xF-xU)/(xD-xU);
	PhiAC=(PhiC-PhiU)/(PhiD-PhiU);
	if(PhiD==PhiU){
		PhiAC=0;
	}
	PhiAF=xAF+(PhiAC-xAC)*(xAF*(xAF-1))/(xAC*(xAC-1));
	PhinHRS=PhiAF*(PhiD-PhiU)+PhiU;
	if(i==0){
		PhiBoundary(PhiwHRS);
	}
	if(i==N-1){
		PhiBoundary(PhieHRS);
	}
	if(j==M-1){
		PhiBoundary(PhinHRS);
	}
	if(j==0&&i<=(N-1)/2){
		PhiInlet(x, i, PhisHRS);
	}
	if(j==0&&i>(N-1)/2){
		PhisHRS=Phi[i][j];
	}
}

void CoefficientbP1(int i, int j, array2D mw, array2D me, array2D ms, array2D mn, double PhiwUDS, double PhieUDS, double PhisUDS, double PhinUDS, double PhiwHRS, double PhieHRS, double PhisHRS, double PhinHRS, array2D& bP){
	bP[i][j]=mw[i][j]*(PhiwHRS-PhiwUDS)-me[i][j]*(PhieHRS-PhieUDS)+ms[i][j]*(PhisHRS-PhisUDS)-mn[i][j]*(PhinHRS-PhinUDS);
}

void CoefficientbP2(int i, int j, array2D rho, array2D V, array2D Phi0, double dt, double Sp, array2D& bP){
	bP[i][j]=bP[i][j]+rho[i][j]*V[i][j]*Phi0[i][j]/dt+Sp*V[i][j];
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

