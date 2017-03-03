//#include <iostream>
#include <fstream>
using namespace std;

//--------------------
//Definici� de constants i variables
//--------------------
//Dades geom�triques
const double B=1; //Profunditat
const double L=1.10; //Dimsnsi� horitzontal
const double H=0.80; //Dimensi� vertical
const double P1x=0.50; //Coordenada x del punt P1
const double P1y=0.40; //Coordenada y del punt P1
const double P2x=0.50; //Coordenada x del punt P2
const double P2y=0.70; //Coordenada y del punt P2
const double L1=P1x; //Longitud del primer tram horitzontal
const double L2=L-P1x; //Longitud del segon tram horitzontal
const double H1=P1y; //Longitud del primer tram vertical
const double H2=P2y-P1y; //Longitud del segon tram vertical
const double H3=H-P2y; //Longitud del tercer tram vertical

//Dades del material
const double rho1=1500.00; //Densitat del material 1
const double rho2=1600.00; //Densitat del material 2
const double rho3=1900.00; //Densitat del material 3
const double rho4=2500.00; //Densitat del material 4
const double cp1=750.00; //Calor espec�fic del material 1
const double cp2=770.00; //Calor espec�fic del material 2
const double cp3=810.00; //Calor espec�fic del material 3
const double cp4=930.00; //Calor espec�fic del material 4
const double lambda1=170.00; //Coeficient de transfer�ncia de calor del material 1
const double lambda2=140.00; //Coeficient de transfer�ncia de calor del material 2
const double lambda3=200.00; //Coeficient de transfer�ncia de calor del material 3
const double lambda4=140.00; //Coeficient de transfer�ncia de calor del material 4

//Dades de contorn
const double TBot=23.00; //Temperatura de la cara inferior
const double qTop=60.00; //Flux de calor per la cara superior
const double TExt=33.00; //Temperatura de l'aire ambient
const double alphaExt=9.00; //Coeficient de transfer�ncia de calor per convecci� de l'aire ambient
const double T0Right=8.00; //Temperatura base de la cara dreta
const double kRight=0.005; //Coeficient amb que var�a la temperatura de la cara dreta amb el temps
const double qv=0.00; //Calor generat per unitat de volum

//Dades num�riques
const int N1=50; //Nombre d'elements en la direcci� horitzontal del primer segment horitzontal
const int N2=60; //Nombre d'elements en la direcci� horitzontal del segon segment horitzontal
const int M1=40; //Nombre de divisions en la direcci� vertical del primer segment vertical
const int M2=30; //Nombre de divisions en la direcci� vertical del segon segment vertical
const int M3=10; //Nombre de divisions en la direcci� vertical del tercer segment vertical
const int Time=10000; //Nombre de divisions en el temps
const double TimeEnd=10000; //Temps final del c�lcul
const double TInit=8.00; //Temperatura inicial del cos
const double delta=0.000001; //Marge d'error per a les iteracions
const double beta=0.5; //Coeficient per determinar el tipus d'integraci� (impl�tica,expl�cuta,intermedia)
const double fr=1.20; //Factor de relaxaci�

//Dades dels punts per exportar
const int output=1; //Valor que indica si s'exporten o no les dades (1 si es volen exportar, 0 si no es volen exportar)
const double xT1=0.65; //Posici� x del primer punt
const double yT1=0.56; //Posici� x del segon punt
const double xT2=0.74; //Posici� y del primer punt
const double yT2=0.72; //Posici� y del segon punt

//C�lculs previs
const double dx1=L1/double(N1); //Dist�ncia horitzontal entre nodes del primer segment horitzontal
const double dx2=L2/double(N2); //Dist�ncia horitzontal entre nodes del segon segment horitzontal
const double dy1=H1/double(M1); //Dist�ncia vertical entre nodes del primer segment vertical
const double dy2=H2/double(M2); //Dist�ncia vertical entre nodes del segon segment vertical
const double dy3=H3/double(M3); //Dist�ncia vertical entre nodes del tercer segment vertical
const double dt=TimeEnd/double(Time); //Interval de temps

//Definici� de variables
double x[N1+N2]; //Definici� del vector de posicions x
double y[M1+M2+M3]; //Definici� del vector de posicions y
double dx[N1+N2]; //Definici� del vector de dist�ncia horitzontal entre nodes
double dy[M1+M2+M3]; //Definici� del vector de dist�nciavertical entre nodes
double Sx[M1+M2+M3]; //Definici� de superf�cies laterals
double Sy[N1+N2]; //Definici� del vector de superf�cies inferior i superior
double V[N1+N2][M1+M2+M3]; //Definici� de la matriu de volums
double rho[N1+N2][M1+M2+M3]; //Definici� de la matriu de la densitat de cada node
double cp[N1+N2][M1+M2+M3]; //Definici� de la matriu del calor espec�fic de cada node
double lambda[N1+N2][M1+M2+M3]; //Definici� de la matriu del coeficient de transfer�ncia de calor de cada node
double T[N1+N2][M1+M2+M3][Time+1]; //Definici� de la matriu de la temperatura de cada node
double lambdaW,lambdaE,lambdaS,lambdaN; //Definici� de les lambdes usades en cada volum finit
double dPW,dPE,dPS,dPN; //Definici� de les dist�ncies entre nodes
double aW,aE,aS,aN,aP,bP; //Definici� de les constants de l'equaci� de la temperatura en un node
double bP1,bP2,bP3,bP4,bP5; //Definici� de constants auxiliars per calcular bP
double time; //Definici� del temps
double maxim; //Definici� del m�xim de la difer�ncia entre la temperatura calculada i la estimada
double TEst[N1+N2][M1+M2+M3]; //Definici� de la temperatura estimada
double T1,T2; //Definici� de la temperatura dels nodes indicats
int iT1, jT1, iT2, jT2; //Definici� dels �ndex de la temperatura del node per exportar
double Titer1, Titer2; //Definici� dels punts intermedis per interpolar

//--------------------
//Funci� main
//--------------------
int main(){
	cout<<"PROGRAMA INICIAT"<<endl<<endl;
	cout<<"Nombre d'elements en la direccio horitzontal: "<<N1+N2<<endl;
	cout<<"Nombre d'elements en la direccio vertical: "<<M1+M2+M3<<endl;
	cout<<"Temps final: "<<TimeEnd<<" s"<<endl;
	cout<<"Intervals de temps: "<<dt<<" s"<<endl;
	cout<<"Factor de relaxacio: "<<fr<<endl;
	if(beta==0){
		cout<<"Esquema explicit"<<endl;
	}
	else if(beta==1){
		cout<<"Esquema implitic"<<endl;
	}
	else if(beta==0.5){
		cout<<"Esquema Crank-Nickolson"<<endl;
	}
	else{
		cout<<"Esquema amb beta="<<beta<<endl;
	}
	cout<<endl;
	
	//Assignaci� de la posici� dels nodes i de les superf�cies de cada eix
	cout<<"Inicialitzant variables..."<<endl;
	for(int i=0;i<=N1-1;i++){
		x[i]=dx1/2+dx1*i; //Coordenada x del node
		dx[i]=dx1; //Dimensi� horitzontal del node
		Sy[i]=dx1*B; //Superf�cie superior i inferior del node
	}
	for(int i=N1;i<=N1+N2-1;i++){
		x[i]=L1+dx2/2+dx2*(i-N1); //Coordenada x del node
		dx[i]=dx2;  //Dimensi� horitzontal del node
		Sy[i]=dx2*B; //Superf�cie superior i inferior del node
	}
	for(int j=0;j<=M1-1;j++){
		y[j]=dy1/2+dy1*j; //Coordenada x del node
		dy[j]=dy1; //Dimensi� vertical del node
		Sx[j]=dy1*B; //Superf�cie lateral del node
	}
	for(int j=M1;j<=M1+M2-1;j++){
		y[j]=H1+dy2/2+dy2*(j-M1); //Coordenada x del node
		dy[j]=dy2; //Dimensi� vertical del node
		Sx[j]=dy2*B; //Superf�cie lateral del node
	}
	for(int j=M1+M2;j<=M1+M2+M3-1;j++){
		y[j]=H1+H2+dy3/2+dy3*(j-M1-M2); //Coordenada x del node
		dy[j]=dy3; //Dimensi� vertical del node
		Sx[j]=dy3*B; //Superf�cie lateral del node
	}
	
	//Assignaci� dels volums de cada volum finit
	for(int i=0;i<=N1+N2-1;i++){
		for(int j=0;j<=M1+M2+M3-1;j++){
			if(x[i]<=L1&&y[j]<=H1){
				V[i][j]=dx1*dy1*B; //Volum del volum finit
			}
			else if(x[i]>L1&&y[j]<=H1){
				V[i][j]=dx2*dy1*B; //Volum del volum finit
			}
			else if(x[i]<=L1&&y[j]>H1&&y[j]<=H1+H2){
				V[i][j]=dx1*dy2*B; //Volum del volum finit
			}
			else if(x[i]>L1&&y[j]>H1&&y[j]<=H1+H2){
				V[i][j]=dx2*dy2*B; //Volum del volum finit
			}
			else if(x[i]<=L1&&y[j]>H1+H2){
				V[i][j]=dx1*dy3*B; //Volum del volum finit
			}
			else if(x[i]>L1&&y[j]>H1+H2){
				V[i][j]=dx2*dy3*B; //Volum del volum finit
			}
			else{
				cout<<"Error. Assignant volum a un punt no contemplat"<<endl;
			}
		}
	}
	
	//Assignaci� de les propietats a cada node
	for(int i=0;i<=N1+N2-1;i++){
		for(int j=0;j<=M1+M2+M3-1;j++){
			if(x[i]<=L1&&y[j]<=H1){ //Nodes que formen part del Material 1
				rho[i][j]=rho1; //Assignaci� de la densitat del material 1 als nodes que formin part del material 1
				cp[i][j]=cp1; //Assignaci� de la densitat del material 1 als nodes que formin part del material 1
				lambda[i][j]=lambda1; //Assignaci� de la densitat del material 1 als nodes que formin part del material 1
			}
			else if(x[i]>L1&&y[j]<=H1+H2){ //Nodes que formen part del Material 2
				rho[i][j]=rho2; //Assignaci� de la densitat del material 2 als nodes que formin part del material 2
				cp[i][j]=cp2; //Assignaci� de la densitat del material 2 als nodes que formin part del material 2
				lambda[i][j]=lambda2; //Assignaci� de la densitat del material 2 als nodes que formin part del material 2
			}
			else if(x[i]<=L1&&y[j]>H1){ //Nodes que formen part del Material 3
				rho[i][j]=rho3; //Assignaci� de la densitat del material 3 als nodes que formin part del material 3
				cp[i][j]=cp3; //Assignaci� de la densitat del material 3 als nodes que formin part del material 3
				lambda[i][j]=lambda3; //Assignaci� de la densitat del material 3 als nodes que formin part del material 3
			}
			else if(x[i]>L1&&y[j]>H1+H2){ //Nodes que formen part del Material 4
				rho[i][j]=rho4; //Assignaci� de la densitat del material 4 als nodes que formin part del material 4
				cp[i][j]=cp4; //Assignaci� de la densitat del material 4 als nodes que formin part del material 4
				lambda[i][j]=lambda4; //Assignaci� de la densitat del material 4 als nodes que formin part del material 4
			}
			else{
				cout<<"Error. Assignant propietats a un punt no contemplat"<<endl;
			}
		}
	}
	
	//Assignaci� de la temperatura inicial a cada node
	for(int i=0;i<=N1+N2-1;i++){
		for(int j=0;j<=M1+M2+M3-1;j++){
			T[i][j][0]=TInit; //Temperatura inicial en el temps 0
			T[i][j][1]=TInit; //Temperatura inicial en el temps 0+dt
			TEst[i][j]=TInit; //Temperatura estimada
		}
	}
	cout<<"Fet"<<endl<<endl;
	
	//C�lcul de la temperatura
	cout<<"Calculant temperatures..."<<endl;
	time=0; //Inicialitzaci� del temps
	for(int t=0;t<=Time-1;t++){
		maxim=1; //Inicialitzaci� del m�xim de la difer�ncia entre la temperatura calculada i la estimada
		if(t==1*Time/10){
			cout<<"Progres: 10%"<<endl;
		}
		if(t==2*Time/10){
			cout<<"Progres: 20%"<<endl;
		}
		if(t==3*Time/10){
			cout<<"Progres: 30%"<<endl;
		}
		if(t==4*Time/10){
			cout<<"Progres: 40%"<<endl;
		}
		if(t==5*Time/10){
			cout<<"Progres: 50%"<<endl;
		}
		if(t==6*Time/10){
			cout<<"Progres: 60%"<<endl;
		}
		if(t==7*Time/10){
			cout<<"Progres: 70%"<<endl;
		}
		if(t==8*Time/10){
			cout<<"Progres: 80%"<<endl;
		}
		if(t==9*Time/10){
			cout<<"Progres: 90%"<<endl;
		}
		while(maxim>delta){
			maxim=0; //M�xim a 0
			for(int i=0;i<=N1+N2-1;i++){
				for(int j=0;j<=M1+M2+M3-1;j++){
					if(i>0&&i<N1+N2-1&&j>0&&j<M1+M2+M3-1){ //Nodes interiors
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aW+aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP4+bP5; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*T[i+1][j][t+1]+aS*T[i][j-1][t+1]+aN*T[i][j+1][t+1]+bP)/aP; //Constant independent
					}
					else if(i==0&&j>0&&j<M1+M2+M3-1){ //Nodes de l'esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt+beta*alphaExt*Sx[j]; //Constant de la temperatura del node actual
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]+(1-beta)*alphaExt*Sx[j]*(TExt-T[i][j][t]); //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+beta*alphaExt*Sx[j]*TExt+bP2+bP3+bP4+bP5; //Constant independent
						T[i][j][t+1]=(aE*T[i+1][j][t+1]+aS*T[i][j-1][t+1]+aN*T[i][j+1][t+1]+bP)/aP; //Equaci� de la temperatura del node
					}
					else if(i==N1+N2-1&&j>0&&j<M1+M2+M3-1){ //Nodes de la dreta
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2; //Dist�ncia entre el node actual i l'extrem dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=lambda[i][j]; //Lambda del node actual
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aW+aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*((T0Right+kRight*time)-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP4+bP5; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*(T0Right+kRight*(time+dt))+aS*T[i][j-1][t+1]+aN*T[i][j+1][t+1]+bP)/aP; //Constant independent
					}
					else if(j==0&&i>0&&i<N1+N2-1){ //Nodes inferiors
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2; //Dist�ncia entre el node actual i l'extrem inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=lambda[i][j]; //Lambda del node actual
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aW+aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(TBot-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP4+bP5; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*T[i+1][j][t+1]+aS*TBot+aN*T[i][j+1][t+1]+bP)/aP; //Constant independent
					}
					else if(j==M1+M2+M3-1&&i>0&&i<N1+N2-1){ //Nodes superiors
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aP=aW+aE+aS+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP5+beta*qTop*dx[i]+(1-beta)*qTop*dx[i]; //Constant independent
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*T[i+1][j][t+1]+aS*T[i][j-1][t+1]+bP)/aP; //Equaci� de la temperatura del node
					}
					else if(i==0&&j==0){ //Node inferior esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2; //Dist�ncia entre el node actual i l'extrem inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=lambda[i][j]; //Lambda del node actual
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt+beta*alphaExt*Sx[j]; //Constant de la temperatura del node actual
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(TBot-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]+(1-beta)*alphaExt*Sx[j]*(TExt-T[i][j][t]); //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+beta*alphaExt*Sx[j]*TExt+bP2+bP3+bP4+bP5; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aE*T[i+1][j][t+1]+aS*TBot+aN*T[i][j+1][t+1]+bP)/aP; //Constant independent
					}
					else if(i==0&&j==M1+M2+M3-1){ //Node superior esquerra
						dPE=dx[i]/2+dx[i+1]/2; //Dist�ncia entre el node actual i el node dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						lambdaE=dPE/(0.5*dx[i]/lambda[i][j]+0.5*dx[i+1]/lambda[i+1][j]); //Mitja harm�nica de la lambda del node actual i la del node dret
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aP=aE+aS+rho[i][j]*V[i][j]*cp[i][j]/dt+beta*alphaExt*Sx[j]; //Constant de la temperatura del node actual
						bP2=(1-beta)*(T[i+1][j][t]-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]+(1-beta)*alphaExt*Sx[j]*(TExt-T[i][j][t]); //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+beta*alphaExt*Sx[j]*TExt+bP2+bP3+bP5+beta*qTop*dx[i]+(1-beta)*qTop*dx[i]; //Constant inependent
						T[i][j][t+1]=(aE*T[i+1][j][t+1]+aS*T[i][j-1][t+1]+bP)/aP; //Equaci� de la temperatura del node
					}
					else if(i==N1+N2-1&&j==0){ //Node inferior dret
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2; //Dist�ncia entre el node actual i l'extrem dret
						dPS=dy[j]/2; //Dist�ncia entre el node actual i l'extrem inferior
						dPN=dy[j]/2+dy[j+1]/2; //Dist�ncia entre el node actual i el node superior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=lambda[i][j]; //Lambda del node actual
						lambdaS=lambda[i][j]; //Lambda del node actual
						lambdaN=dPN/(0.5*dy[j]/lambda[i][j]+0.5*dy[j+1]/lambda[i][j+1]); //Mitja harm�nica de la lambda del node actual i la del node superior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aN=beta*lambdaN*Sy[i]/dPN; //Constant de la temperatura del node superior
						aP=aW+aE+aS+aN+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*((T0Right+kRight*time)-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(TBot-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP4=(1-beta)*(T[i][j+1][t]-T[i][j][t])*lambdaN*Sy[i]/dPN; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP4+bP5; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*(T0Right+kRight*(time+dt))+aS*TBot+aN*T[i][j+1][t+1]+bP)/aP; //Constant independent
					}
					else if(i==N1+N2-1&&j==M1+M2+M3-1){ //Node superior dret
						dPW=dx[i]/2+dx[i-1]/2; //Dist�ncia entre el node actual i el node esquerra
						dPE=dx[i]/2; //Dist�ncia entre el node actual i l'extrem dret
						dPS=dy[j]/2+dy[j-1]/2; //Dist�ncia entre el node actual i el node inferior
						lambdaW=dPW/(0.5*dx[i]/lambda[i][j]+0.5*dx[i-1]/lambda[i-1][j]); //Mitja harm�nica de la lambda del node actual i la del node esquerra
						lambdaE=lambda[i][j]; //Lambda del node actual
						lambdaS=dPS/(0.5*dy[j]/lambda[i][j]+0.5*dy[j-1]/lambda[i][j-1]); //Mitja harm�nica de la lambda del node actual i la del node inferior
						aW=beta*lambdaW*Sx[j]/dPW; //Constant de la temperatura del node esquerra
						aE=beta*lambdaE*Sx[j]/dPE; //Constant de la temperatura del node dret
						aS=beta*lambdaS*Sy[i]/dPS; //Constant de la temperatura del node inferior
						aP=aW+aE+aS+rho[i][j]*V[i][j]*cp[i][j]/dt; //Constant de la temperatura del node actual
						bP1=(1-beta)*(T[i-1][j][t]-T[i][j][t])*lambdaW*Sx[j]/dPW; //Constant auxiliar per calcular bP
						bP2=(1-beta)*((T0Right+kRight*time)-T[i][j][t])*lambdaE*Sx[j]/dPE; //Constant auxiliar per calcular bP
						bP3=(1-beta)*(T[i][j-1][t]-T[i][j][t])*lambdaS*Sy[i]/dPS; //Constant auxiliar per calcular bP
						bP5=beta*qv*V[i][j]+(1-beta)*qv*V[i][j]; //Constant auxiliar per calcular bP
						bP=rho[i][j]*V[i][j]*cp[i][j]*T[i][j][t]/dt+bP1+bP2+bP3+bP5+beta*qTop*dx[i]+(1-beta)*qTop*dx[i]; //Equaci� de la temperatura del node
						T[i][j][t+1]=(aW*T[i-1][j][t+1]+aE*(T0Right+kRight*(time+dt))+aS*T[i][j-1][t+1]+bP)/aP; //Constant independent
					}
					else{
						cout<<"Error. Calculant temperatura a un punt no contemplat"<<endl;
					}
					
					//Comprovaci� de la converg�ncia
					if(T[i][j][t+1]-TEst[i][j]>=0){
						if(T[i][j][t+1]-TEst[i][j]>maxim){
							maxim=T[i][j][t+1]-TEst[i][j]; //Assignaci� al m�xim el valor de la difer�ncia entre la temperatura estimada i la calculada
						}
					}
					else{
						if(TEst[i][j]-T[i][j][t+1]>maxim){
							maxim=TEst[i][j]-T[i][j][t+1]; //Assignaci� al m�xim el valor de la difer�ncia entre la temperatura estimada i la calculada
						}
					}
					TEst[i][j]=TEst[i][j]+fr*(T[i][j][t+1]-TEst[i][j]); //La temperatura calculada passa a ser la temperatura estimada per a la seg�ent iteraci�
				}
			}
		} //Fi del bucle iteratiu
		time=time+dt; //Seg�ent increment de temps
	}
	cout<<"Fet"<<endl<<endl;
	
	//C�lcul de la posici� dels punts d'on s'exportar� la temperatura
	for(int i=0;i<=N1+N2-2;i++){
		if(x[i]<=xT1&&x[i+1]>xT1){
			iT1=i; //Assignaci� de la posici� del primer punt el el vector x
		}
		if(x[i]<=xT2&&x[i+1]>xT2){
			iT2=i; //Assignaci� de la posici� del segon punt el el vector x
		}
	}
	for(int j=0;j<=M1+M2+M3-2;j++){
		if(y[j]<=yT1&&y[j+1]>yT1){
			jT1=j; //Assignaci� de la posici� del primer punt el el vector y
		}
		if(y[j]<=yT2&&y[j+1]>yT2){
			jT2=j; //Assignaci� de la posici� del segon punt el el vector y
		}
	}
	
	//Exportar temperatures a un document txt
	if(output==1){
		cout<<"Exportant dades..."<<endl;
		ofstream file("SERRA.dat"); //Creaci� del fitxer .txt on guardar les temperatures
		time=0; //Inicialitzaci� del temps
		for(int t=0;t<=Time;t++){
			Titer1=T[iT1][jT1][t]+(T[iT1+1][jT1][t]-T[iT1][jT1][t])*(xT1-x[iT1])/(x[iT1+1]-x[iT1]); //Interpolaci� en x per j per a la temperatura del primer punt
			Titer2=T[iT1][jT1+1][t]+(T[iT1+1][jT1+1][t]-T[iT1][jT1+1][t])*(xT1-x[iT1])/(x[iT1+1]-x[iT1]); //Interpolaci� en x per j+1 per a la temperatura del primer punt
			T1=Titer1+(Titer2-Titer1)*(yT1-y[jT1])/(y[jT1+1]-y[jT1]); //Interpolaci� final de la temperatura del primer punt
			Titer1=T[iT2][jT2][t]+(T[iT2+1][jT2][t]-T[iT2][jT2][t])*(xT2-x[iT2])/(x[iT2+1]-x[iT2]); //Interpolaci� en x per j per a la temperatura del segon punt
			Titer2=T[iT2][jT2+1][t]+(T[iT2+1][jT2+1][t]-T[iT2][jT2+1][t])*(xT2-x[iT2])/(x[iT2+1]-x[iT2]); //Interpolaci� en x per j+1 per a la temperatura del segon punt
			T2=Titer1+(Titer2-Titer1)*(yT2-y[jT2])/(y[jT2+1]-y[jT2]); //Interpolaci� final de la temperatura del primer punt
			file<<time<<"	"<<T1<<"	"<<T2<<endl; //Valors de temps i temperatura delsnodes indicats
			time++; //Seg�ent instant
		}
		file.close();
		cout<<"Fet"<<endl<<endl;
	}
	
	cout<<"PROGRAMA ACABAT"<<endl;
	
	return 0;
}
