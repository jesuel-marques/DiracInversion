#include <stdio.h>
#include <math.h>
#include "../ParametrosSU3.h"
#include "GlobalSU3.h"

#include "FuncoesMatematicas.h"

double identidadeSU3[3][3][2];	//	A matriz identidade 3x3. Inicializada em InicializarUSU3 no arquivo RedeSU3.c
double nuloSU3[3][3][2];	//	A matriz 0 complexa 3x3. Inicializada em InicializarUSU3 no arquivo RedeSU3.c

// double lambda[9][3][3][2];	//	Os geradores da álgebra de SU(3). Inicializados em InicializarLambda nesse arquivo.

void CopiarSU3(double u[3][3][2], double ucopiado[3][3][2]){

// Copia matriz u para matriz ucopiado

	for(int a=0;a<3;a++)
		for(int b=0;b<3;b++)
			CopiarComplexo(u[a][b],ucopiado[a][b]);
}

void SomaSU3(double u[3][3][2],double resultado[3][3][2]){

//	Soma u com resultado e acumula em resultado

	for(int a=0;a<3;a++)
		for(int b=0;b<3;b++)
			SomaComplexo(u[a][b],resultado[a][b]);
}

void DiferencaSU3(double u[3][3][2], double v[3][3][2], double umenosv[3][3][2]){

//	Calcula a diferença entre a matriz u e a matriz v e coloca resultado em umenosv

	for(int a=0;a<3;a++)
		for(int b=0;b<3;b++)
			DiferencaComplexo(u[a][b],v[a][b],umenosv[a][b]);
		
}

void TrSU3(double u[3][3][2],double tr[2]){

//	Calcula o traço da matriz u e coloca resultado, um número complexo, em tr

	CopiarComplexo(ZeroComplexo,tr);
	for(int a=0;a<3;a++)
		SomaComplexo(u[a][a],tr);
}

double ReTr(double u[3][3][2]){

//	Devolve a parte real do traço da matriz u

	double tr[2];
	TrSU3(u,tr);
	return tr[0];
}

void DeterminanteSU3(double u[3][3][2], double det[2]){

//	Calcula o determinante da matriz u

	double aux1[2], aux2[2], aux3[2], aux4[2];

	ProdutoComplexo(u[1][1],u[2][2],aux1);
	ProdutoComplexo(u[1][2],u[2][1],aux2);
	DiferencaComplexo(aux1,aux2,aux3);
	ProdutoComplexo(u[0][0],aux3,det);

	ProdutoComplexo(u[1][2],u[2][0],aux1);
	ProdutoComplexo(u[1][0],u[2][2],aux2);
	DiferencaComplexo(aux1,aux2,aux3);
	ProdutoComplexo(u[0][1],aux3,aux4);
	SomaComplexo(aux4,det);

	ProdutoComplexo(u[1][0],u[2][1],aux1);
	ProdutoComplexo(u[1][1],u[2][0],aux2);
	DiferencaComplexo(aux1,aux2,aux3);
	ProdutoComplexo(u[0][2],aux3,aux4);
	SomaComplexo(aux4,det);
}

void ConjHermSU3(double X[3][3][2],double Xch[3][3][2]){

// Calcula conjugado hermitiano da matriz X e coloca em Xch

	for(int i=0;i<3;i++)
		for(int j=0;j<=i;j++){

				ConjugadoComplexo(X[i][j],Xch[j][i]);
				ConjugadoComplexo(X[j][i],Xch[i][j]);
	
		}
}

void ProdutoSU3(double Y[3][3][2], double Z[3][3][2], double YZ[3][3][2]){

// Calcula produto de 2 matrizes SU(3) Y e Z e coloca resultado em YZ

	CopiarSU3(nuloSU3,YZ);

	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		for(int k=0;k<3;k++){	// índice mudo
			YZ[i][j][RE]+=Y[i][k][RE]*Z[k][j][RE]-Y[i][k][IM]*Z[k][j][IM];
			YZ[i][j][IM]+=Y[i][k][RE]*Z[k][j][IM]+Y[i][k][IM]*Z[k][j][RE];
		}
}

void ProdutoSU3Tres(double X[3][3][2], double Y[3][3][2], double Z[3][3][2], double XYZ[3][3][2]){

// Calcula produto de 3 matrizes SU(3) X, Y e Z e coloca resultado em XYZ

	double XY[3][3][2];
	ProdutoSU3(X,Y,XY);
	ProdutoSU3(XY,Z,XYZ);
}

void ProdutoSU3Quatro(double W[3][3][2], double X[3][3][2], double Y[3][3][2], double Z[3][3][2], double WXYZ[3][3][2]){

// Calcula produto de 4 matrizes SU(3) W, X, Y e Z e coloca resultado em WXYZ

	double WXY[3][3][2];
	ProdutoSU3Tres(W,X,Y,WXY);
	ProdutoSU3(WXY,Z,WXYZ);
}

void PotenciaSU3(double u[3][3][2],double unapotencias[kMAXbinomial+1][3][3][2]){

//	Calcula as potências de u e as coloca em diferentes linhas de unapotencias.

	for(int a=0;a<=3;a++)
	for(int b=0;b<=3;b++){
		CopiarComplexo(unapotencias[0][a][b],identidadeSU3[a][b]); // u elevado a 0 é a identidade
		CopiarComplexo(unapotencias[1][a][b],u[a][b]);	// u elevadoa 1 é si mesmo
	}

	for(int i=2;i<=kMAXbinomial;i++)
		ProdutoSU3(unapotencias[i-1],u,unapotencias[i]);
	
}

void MultiplicacaoEscalarSU3(double u[3][3][2],double num[2], double numvezesu[3][3][2]){

//	Função que retorna a matriz SU(3) multiplicada por um escalar complexo.

	for(int a=0;a<3;a++)
	for(int b=0;b<3;b++)
		ProdutoComplexo(num,u[a][b],numvezesu[a][b]);	//	Basta multiplicar cada componente.

}

void ProjecaoSU3(double u[3][3][2], double uSU3[3][3][2]){	

//	Projeta matriz u para o grupo SU(3) produzindo a matriz uSU3, de determinante 1 e cuja inversa é a tranposta conjugada. 
//	Segue método explicado em Gattringer ao redor da Eq. 4.27

	double somamodulo=0.0;
 	double umsobremodulo[2]={0.0,0.0};
	double unovoestrela[3][2];
	double vnovoestrela[3][2];
	double termovunovoestrela[2]={0.0,0.0};
	double vunovoestrela[2]={0.0,0.0};
	double unovovunovoestrela[3][2];
	double vlinha[3][2];

	double auxprod1[2];
	double auxprod2[2];

	for(int i=0;i<3;i++)
		somamodulo+=ModuloQuadComplexo(u[0][i]);

	umsobremodulo[RE]=1.0/sqrt(somamodulo);

	for(int i=0;i<3;i++){
		ProdutoComplexo(umsobremodulo,u[0][i],uSU3[0][i]);	//	Calcula unovo, primeira linha da matriz projetada a partir de u simplesmente dividindo pelo módulo
		ConjugadoComplexo(uSU3[0][i],unovoestrela[i]);	//	unovoestrela é o conjugado de unovo
		ProdutoComplexo(u[1][i],unovoestrela[i],termovunovoestrela);	//	Termo da projeção de v sobre unovoestrela
		SomaComplexo(termovunovoestrela,vunovoestrela);
	}

	for(int i=0;i<3;i++){	//	Método de Gram-Schidmt
		ProdutoComplexo(uSU3[0][i],vunovoestrela,unovovunovoestrela[i]);	//	termo que será subtraído de v no método de Gram-Schmidt
		DiferencaComplexo(u[1][i],unovovunovoestrela[i],vlinha[i]);
	}


	somamodulo=0.0;

	for(int i=0;i<3;i++){
		somamodulo+=ModuloQuadComplexo(vlinha[i]);	//	Segunda linha da matriz projetada antes da projeção
	}

	umsobremodulo[RE]=1.0/sqrt(somamodulo);	

	for(int i=0;i<3;i++){
		ProdutoComplexo(umsobremodulo,vlinha[i],uSU3[1][i]);	//	Calcula vnovo, segunda linha da matriz projetada a partir de vlinha
		ConjugadoComplexo(uSU3[1][i],vnovoestrela[i]);	//	vnovoestrela é o conjugado de vnovo
	}

	//	A terceira linha é calculada tirando o produto cruzado de unovoestrela com vnovoestrela

	ProdutoComplexo(unovoestrela[1],vnovoestrela[2],auxprod1);
	ProdutoComplexo(unovoestrela[2],vnovoestrela[1],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,uSU3[2][0]);	//	primeira componente

	
	ProdutoComplexo(unovoestrela[2],vnovoestrela[0],auxprod1);
	ProdutoComplexo(unovoestrela[0],vnovoestrela[2],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,uSU3[2][1]);	//	segunda componente


	ProdutoComplexo(unovoestrela[0],vnovoestrela[1],auxprod1);
	ProdutoComplexo(unovoestrela[1],vnovoestrela[0],auxprod2);
	DiferencaComplexo(auxprod1,auxprod2,uSU3[2][2]);	//	terceira componente

}

void DecomporAlgebraSU3(double g[3][3][2],double ga[9]){ //	Decompõe nas componentes da álgebra de SU(3) usando convenções do Gattringer.

	// Fórmulas para coeficientes obtidas no Mathematica.
	
	ga[1]=(g[0][1][RE]+g[1][0][RE]);
	
	ga[2]=(-g[0][1][IM]+g[1][0][IM]);
	
	
	ga[3]=(g[0][0][RE]-g[1][1][RE]);

	
	ga[4]=(g[0][2][RE]+g[2][0][RE]);
	
	ga[5]=(-g[0][2][IM]+g[2][0][IM]);

	
	ga[6]=(g[1][2][RE]+g[2][1][RE]);
	
	ga[7]=(-g[1][2][IM]+g[2][1][IM]);

	
	ga[8]=(g[0][0][RE]+g[1][1][RE]-2.0*g[2][2][RE])/(pow(3.0,0.5));


}

void ImprimirMatrizSU3(double M[3][3][2],char *nome ){

	//	Imprime matriz M com formato para copiar e colar no Mathematica
	
	printf("\n\n %s \n",nome);

		printf("\n");


		printf("{");
		for(int i=0;i<3;i++){
			printf("{");

		for(int j=0;j<3;j++){
			printf("%.16lf+I(%.16lf)",M[i][j][RE],M[i][j][IM]);
			if(j!=2) printf(",");
		}
		if(i!=2){

			printf("},\n");
		}	
		}

		printf("}}");

		getchar();

}