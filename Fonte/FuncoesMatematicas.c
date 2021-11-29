#include <math.h>

#include "GlobalSU3.h"
#include "../ParametrosSU3.h"

#include "SU3.h"

#include <stdbool.h>

double ZeroComplexo[2]={0.0,0.0};
double UmComplexo[2]={1.0,0.0};

//Funções Matemáticas

double Sign(double  x){

//	Função que retorna o sinal de um número x (+1 se positivo e -1 se negativo)

	return ((x > 0) - (x < 0)); 

}

double pow2(double x){ 

//	Função que retorna o quadrado de um número

	return x*x;

}

double  ProdutoInterno(double  x[4],double  u[4]){ 

//	Função que retorna o produto escalar de dois vetores quadridimensionais representando
// (através da representação de Cayley-Klein) matrizes SU(2). Usada na multiplicação de 
//	dois elementos SU(2).


	double  prodint=0;

	prodint=x[0]*u[0];	//	A componente 0, associada à identidade, tem sinal + no produto.

	for(int b=1;b<=3;b++){
		prodint+=-x[b]*u[b]; // As componentes 1, 2 e 3, associadas às matrizes de Pauli,
							 //	tem sinal - no produto.
	}

	return prodint;

}

void ProdutoExterno(double  x[4], double  u[4], double  prodext[4]){ 

//	Função que retorna o produto vetorial de dois vetores quadridimensionais representando matrizes SU(2).
//	Usada na multiplicação de dois elementos SU(2).
//	De fato o produto é tomado apenas pelas componentes 1, 2 e 3

	prodext[1]=x[2]*u[3]-x[3]*u[2];
	prodext[2]=x[3]*u[1]-x[1]*u[3];
	prodext[3]=x[1]*u[2]-x[2]*u[1];
	
}

void CopiarComplexo(double x[2],double xcopiado[2]){
	xcopiado[RE]=x[RE];
	xcopiado[IM]=x[IM];
}

void ConjugadoComplexo(double x[2],double xestrela[2]){
	xestrela[RE]=x[RE];
	xestrela[IM]=-x[IM];
}

void SomaComplexo(double x[2],double soma[2]){
	soma[RE]+=x[RE];
	soma[IM]+=x[IM];
}

void DiferencaComplexo(double x[2],double y[2],double diferenca[2]){
	diferenca[RE]=x[RE]-y[RE];
	diferenca[IM]=x[IM]-y[IM];
}

void ProdutoComplexo(double x[2], double y[2],double xy[2]){
	xy[RE]=x[RE]*y[RE]-x[IM]*y[IM];
	xy[IM]=x[IM]*y[RE]+x[RE]*y[IM];
}

void ProdutoConjugadoComplexo(double x[2], double y[2],double xestrelay[2]){
	xestrelay[RE]=x[RE]*y[RE]+x[IM]*y[IM];
	xestrelay[IM]=x[RE]*y[IM]-x[IM]*y[RE];
}

void ProdutoComplexoTres(double x[2],double y[2], double w[2],double xyw[2]){
	double xy[2];

	ProdutoComplexo(x,y,xy);
	ProdutoComplexo(xy,w,xyw);
}

void ProdutoComplexoQuatro(double x[2],double y[2], double w[2], double z[2], double xywz[2]){
	double xy[2];
	double wz[2];
	
	ProdutoComplexo(x,y,xy);
	ProdutoComplexo(w,z,wz);
	ProdutoComplexo(xy,wz,xywz);
}

double ModuloQuadComplexo(double x[2]){
	return pow2(x[RE])+pow2(x[IM]);
}

void DivisaoComplexo(double x[2], double y[2], double xpory[2]){
	double yestrela[2], xyestrela[2];
	ConjugadoComplexo(y,yestrela);
	ProdutoComplexo(x,yestrela,xyestrela);
	xpory[RE]=xyestrela[RE]/ModuloQuadComplexo(y);
	xpory[IM]=xyestrela[IM]/ModuloQuadComplexo(y);
}

void Inversa3por3(double a[3][3][2],double ainv[3][3][2]){

	//	Calcula inversa de matriz 3x3 complexa explicitamente

	double auxprod1[2];
	double auxprod2[2];
	double auxsub[2];

	double det[2];

	DeterminanteSU3(a,det);

	ProdutoComplexo(a[1][1],a[2][2],auxprod1);
	ProdutoComplexo(a[1][2],a[2][1],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[0][0]);


	ProdutoComplexo(a[0][2],a[2][1],auxprod1);
	ProdutoComplexo(a[0][1],a[2][2],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[0][1]);


	ProdutoComplexo(a[0][1],a[1][2],auxprod1);
	ProdutoComplexo(a[0][2],a[1][1],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[0][2]);


	ProdutoComplexo(a[1][2],a[2][0],auxprod1);
	ProdutoComplexo(a[1][0],a[2][2],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[1][0]);


	ProdutoComplexo(a[0][0],a[2][2],auxprod1);
	ProdutoComplexo(a[0][2],a[2][0],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[1][1]);


	ProdutoComplexo(a[0][2],a[1][0],auxprod1);
	ProdutoComplexo(a[0][0],a[1][2],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[1][2]);


	ProdutoComplexo(a[1][0],a[2][1],auxprod1);
	ProdutoComplexo(a[1][1],a[2][0],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[2][0]);


	ProdutoComplexo(a[0][1],a[2][0],auxprod1);
	ProdutoComplexo(a[0][0],a[2][1],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[2][1]);


	ProdutoComplexo(a[0][0],a[1][1],auxprod1);
	ProdutoComplexo(a[0][1],a[1][0],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,auxsub);
	DivisaoComplexo(auxsub,det,ainv[2][2]);

}

void Calcularexp(double amomento[d],int posicao[d],double expipx[2]){

	//	Calcula exponencial de (i sum_mu p_mu  x_mu)
	
	double px=0.0;
	for(int i=0;i<d;i++){
		px+=(double)amomento[i]*posicao[i];
	}
	expipx[RE]=cos(px);
	expipx[IM]=-sin(px);

}

void ProdutocomEscalarComplexo4x4(double a[2], double X[4][4][2], double aX[4][4][2]){
	for(int alfa=0;alfa<4;alfa++)
	for(int beta=0;beta<4;beta++)
		ProdutoComplexo(a,X[alfa][beta],aX[alfa][beta]);
}

void DiferencaComplexo4x4(double Y[4][4][2], double Z[4][4][2], double YmenosZ[4][4][2]){
		
	for(int alfa=0;alfa<4;alfa++)
	for(int beta=0;beta<4;beta++)
		DiferencaComplexo(Y[alfa][beta],Z[alfa][beta],YmenosZ[alfa][beta]);
}

void ProdutoComplexo4x4(double Y[4][4][2], double Z[4][4][2], double YZ[4][4][2]){

// Calcula produto de 2 matrizes complexas 4x4 Y e Z e coloca resultado em YZ

	double termo[2];
	
	for(int alfa=0;alfa<4;alfa++)
	for(int beta=0;beta<4;beta++)
		CopiarComplexo(ZeroComplexo,YZ[alfa][beta]);

	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++)
		for(int k=0;k<4;k++){	// índice mudo
			ProdutoComplexo(Y[i][k],Z[k][j],termo);
			SomaComplexo(termo,YZ[i][j]);
		}
}

void ComutadorComplexo4x4(double Y[4][4][2], double Z[4][4][2], double YZmenosZY[4][4][2]){
	double YZ[4][4][2];
	double ZY[4][4][2];

	ProdutoComplexo4x4(Y,Z,YZ);
	ProdutoComplexo4x4(Z,Y,ZY);

	DiferencaComplexo4x4(YZ,ZY,YZmenosZY);
}