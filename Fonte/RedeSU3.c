#include <stdio.h>
#include <time.h>
#include <math.h>

#include "../ParametrosSU3.h"
#include "GlobalSU3.h"
#include "GlobalFixacaoGaugeSU3.h"

#include "FuncoesMatematicas.h"
#include "SU3.h"
#include "ranlux.h"

//	Submatrizes SU(2) de SU(3) para Inicialização Quente, Banho Térmico e Fixação de Calibre por Los Alamos ou Variantes

double R[3][3][2];	
double S[3][3][2];
double T[3][3][2];

void InicializarUSU3(int Partida){	

//Inicializa a rede com uma partida fria ou quente.

	//	Variáveis para Partida Quente.
	double moda, modb;
	double fasea, faseb;

	double auxalea[3];

	double X[3][3][2];

	for(int e=0;e<3;e++)
	for(int f=0;f<3;f++)
		for(int g=RE;g<=IM;g++){

			R[e][f][g]=0.0;
			S[e][f][g]=0.0;
			T[e][f][g]=0.0;
			identidadeSU3[e][f][g]=0.0;
			nuloSU3[e][f][g]=0.0;
		}
	
	
	identidadeSU3[0][0][RE]=1.0;
	identidadeSU3[1][1][RE]=1.0;
	identidadeSU3[2][2][RE]=1.0;

	R[2][2][RE]=1.0;
	S[1][1][RE]=1.0;
	T[0][0][RE]=1.0;

	switch(Partida){
		case FRIA:
			for(int t=0;t<Nt;t++)
				for(int i=0;i<Nxyz;i++)
				for(int j=0;j<Nxyz;j++)
				for(int k=0;k<Nxyz;k++){
					CopiarSU3(identidadeSU3,G[t][i][j][k]);
					for(int mu=0;mu<d;mu++){
						CopiarSU3(identidadeSU3,U[t][i][j][k][mu]);	//	Na partida fria coloca matrizes identidade em U
						CopiarSU3(identidadeSU3,Uaux[t][i][j][k][mu]);	//	Uaux e G são utilizados na fixação de calibre
						
					}
				}
			
			break;

		case QUENTE:

			for(int t=0;t<Nt;t++)
				for(int i=0;i<Nxyz;i++)
				for(int j=0;j<Nxyz;j++)
				for(int k=0;k<Nxyz;k++)
					for(int mu=0;mu<d;mu++){

						ranlxd(auxalea,3);	//	Na partida quente coloca matrizes randomicas pertencentes ao grupo em U
					
						moda=sqrt(auxalea[0]);
						modb=sqrt(1.0-pow2(moda));

						fasea=(2.0*M_PI)*(auxalea[1]);
						faseb=(2.0*M_PI)*(auxalea[2]);

						R[0][0][RE]=moda*cos(fasea);
						R[0][0][IM]=moda*sin(fasea);

						R[0][1][RE]=modb*cos(faseb);
						R[0][1][IM]=modb*sin(faseb);

						R[1][0][RE]=-modb*cos(faseb);
						R[1][0][IM]=modb*sin(faseb);

						R[1][1][RE]=moda*cos(fasea);
						R[1][1][IM]=-moda*sin(fasea);
						
						ranlxd(auxalea,3);
					
						moda=sqrt(auxalea[0]);
						modb=sqrt(1.0-pow2(moda));

						fasea=(2.0*M_PI)*(auxalea[1]);
						faseb=(2.0*M_PI)*(auxalea[2]);
						
						S[0][0][RE]=moda*cos(fasea);
						S[0][0][IM]=moda*sin(fasea);

						S[0][2][RE]=modb*cos(faseb);
						S[0][2][IM]=modb*sin(faseb);

						S[2][0][RE]=-modb*cos(faseb);
						S[2][0][IM]=modb*sin(faseb);

						S[2][2][RE]=moda*cos(fasea);
						S[2][2][IM]=-moda*sin(fasea);

						ranlxd(auxalea,3);
					
						moda=sqrt(auxalea[0]);
						modb=sqrt(1.0-pow2(moda));

						fasea=(2.0*M_PI)*(auxalea[1]);
						faseb=(2.0*M_PI)*(auxalea[2]);

						T[1][1][RE]=moda*cos(fasea);
						T[1][1][IM]=moda*sin(fasea);

						T[1][2][RE]=modb*cos(faseb);
						T[1][2][IM]=modb*sin(faseb);

						T[2][1][RE]=-modb*cos(faseb);
						T[2][1][IM]=modb*sin(faseb);

						T[2][2][RE]=moda*cos(fasea);
						T[2][2][IM]=-moda*sin(fasea);

						ProdutoSU3Tres(R,S,T,X);

						ProjecaoSU3(X,U[t][i][j][k][mu]);
						CopiarSU3(U[t][i][j][k][mu],Uaux[t][i][j][k][mu]);
						CopiarSU3(identidadeSU3,G[t][i][j][k]);

						}
			break;


		default :
			printf("Definir Partida.\n");

	}


}

void TransformacaoCalibreRandomicaSU3(){

	double moda, modb;
	double fasea, faseb;

	double auxalea[3];

	double X[3][3][2];

	for(int t=0;t<Nt;t++)
		for(int i=0;i<Nxyz;i++)
		for(int j=0;j<Nxyz;j++)
		for(int k=0;k<Nxyz;k++){
						
			ranlxd(auxalea,3);	//	Na partida quente coloca matrizes randomicas pertencentes ao grupo em U
					
			moda=sqrt(auxalea[0]);
			modb=sqrt(1.0-pow2(moda));

			fasea=(2.0*M_PI)*(auxalea[1]);
			faseb=(2.0*M_PI)*(auxalea[2]);

			R[0][0][RE]=moda*cos(fasea);
			R[0][0][IM]=moda*sin(fasea);

			R[0][1][RE]=modb*cos(faseb);
			R[0][1][IM]=modb*sin(faseb);

			R[1][0][RE]=-modb*cos(faseb);
			R[1][0][IM]=modb*sin(faseb);

			R[1][1][RE]=moda*cos(fasea);
			R[1][1][IM]=-moda*sin(fasea);
			
			ranlxd(auxalea,3);
		
			moda=sqrt(auxalea[0]);
			modb=sqrt(1.0-pow2(moda));

			fasea=(2.0*M_PI)*(auxalea[1]);
			faseb=(2.0*M_PI)*(auxalea[2]);
			
			S[0][0][RE]=moda*cos(fasea);
			S[0][0][IM]=moda*sin(fasea);

			S[0][2][RE]=modb*cos(faseb);
			S[0][2][IM]=modb*sin(faseb);

			S[2][0][RE]=-modb*cos(faseb);
			S[2][0][IM]=modb*sin(faseb);

			S[2][2][RE]=moda*cos(fasea);
			S[2][2][IM]=-moda*sin(fasea);

			ranlxd(auxalea,3);
		
			moda=sqrt(auxalea[0]);
			modb=sqrt(1.0-pow2(moda));

			fasea=(2.0*M_PI)*(auxalea[1]);
			faseb=(2.0*M_PI)*(auxalea[2]);

			T[1][1][RE]=moda*cos(fasea);
			T[1][1][IM]=moda*sin(fasea);

			T[1][2][RE]=modb*cos(faseb);
			T[1][2][IM]=modb*sin(faseb);

			T[2][1][RE]=-modb*cos(faseb);
			T[2][1][IM]=modb*sin(faseb);

			T[2][2][RE]=moda*cos(fasea);
			T[2][2][IM]=-moda*sin(fasea);

			ProdutoSU3Tres(R,S,T,X);

			ProjecaoSU3(X,G[t][i][j][k]);
		}
}

void CopiarConfigSU3(double ConfigU[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2],double ConfigUcopiado[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2]){
	for(int t=0;t<Nt;t++)
		for(int i=0;i<Nxyz;i++)
		for(int j=0;j<Nxyz;j++)
		for(int k=0;k<Nxyz;k++)
			for(int mu=0; mu<d;mu++){
			CopiarSU3(ConfigU[t][i][j][k][mu],ConfigUcopiado[t][i][j][k][mu]);
		}
}


void CarregarUSU3(char NomeArquivo[ComprimentoMaxArquivo]){	

//	Carrega a partir de um arquivo uma configuração de elos no programa

	FILE *Arquivo;
	
	int posicao[d];
	int mu, a, b, c;
	double elementoleitura;

	InicializarUSU3(FRIA);	//	Antes de carregar coloca matrizes identidade em toda a rede por segurança

	printf("%s\n",NomeArquivo);


	Arquivo=fopen(NomeArquivo,"r");	

	while(fscanf(Arquivo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",&posicao[0],&posicao[1],&posicao[2],&posicao[3],&mu,&a,&b,&c,&elementoleitura)!=EOF){
		U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu][a][b][c]=elementoleitura;
		Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu][a][b][c]=elementoleitura;
	}

	// if(fread(U,sizeof(U),1,Arquivo)==1){
	// 	printf("Config %s carregada.\n",NomeArquivo);
	// }
	// else {
	// 	printf("Erro carregando configuração.");
	// }

	fclose(Arquivo);

		for(int t=0;t<Nt;t++)
				for(int i=0;i<Nxyz;i++)
				for(int j=0;j<Nxyz;j++)
				for(int k=0;k<Nxyz;k++){
					CopiarSU3(identidadeSU3,G[t][i][j][k]);
					for(int mu=0;mu<d;mu++)
						CopiarSU3(U[t][i][j][k][mu],Uaux[t][i][j][k][mu]);
				}
}

void ImprimirUSU3(char NomeArquivo[ComprimentoMaxArquivo]){	

//	Imprime em um arquivo a configuração atual de elos.

	FILE *Arquivo;

	Arquivo=fopen(NomeArquivo,"w+");

	fwrite(U, sizeof(U),1, Arquivo);
	
	// for(int t=0;t<Nt;t++)
	// 	for(int i=0;i<Nxyz;i++)
	// 	for(int j=0;j<Nxyz;j++)
	// 	for(int k=0;k<Nxyz;k++)
	// 		for(int mu=0;mu<d;mu++)
	// 			for(int a=0;a<3;a++)
	// 			for(int b=0;b<3;b++)
	// 				for(int c=RE;c<=IM;c++)
	// 					fprintf(Arquivo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.16lf\n",t,i,j,k,mu,a,b,c,U[t][i][j][k][mu][a][b][c]);

	fclose(Arquivo);
}

void CopiarPosicao(int v[d],int u[d]){
	for(int i=0;i<d;i++)
		u[i]=v[i];
}


void SomaVetoresPosicao(int v[d], int u[d], int vmaisu[d]){

//	Dados dois vetores no espaço real, retorna a sua soma, levando em consideração as condições
//	periódicas de contorno. 	

	if((v[0]+u[0])>=0)	//	Componente temporal
		vmaisu[0]=((v[0]+u[0])%Nt);
	else
		vmaisu[0]=(((v[0]+u[0])%Nt+Nt)%Nt);

	for(int i=1;i<d;i++){ //	Componentes espaciais 

		if((v[i]+u[i])>=0)
			vmaisu[i]=((v[i]+u[i])%Nxyz);
		else
			vmaisu[i]=(((v[i]+u[i])%Nxyz+Nxyz)%Nxyz);

	}	

}

void PosicoesVizinhas(int posicao[d], int mu, int sentidomu, int nu, int sentidonu, int posicaomaismu[d], int posicaomaismumaisnu[d],int posicaomaisnu[d], int posicaomaismumenosnu[d], int posicaomenosnu[d]){ //OK

//	Função que retorna as posições vizinhas, dado um sítio e direções mu e nu.
//	O cálculo é feito utilizando a função que soma a posição dos vetores com os versores nas direções
//	mu e nu.

	int versormu[d]={0,0,0,0};	int versormenosmu[d]={0,0,0,0};
	int versornu[d]={0,0,0,0};	int versormenosnu[d]={0,0,0,0};

	if(sentidomu==FRENTE){
		versormu[mu]=1;	versormenosmu[mu]=-1;
	}
	else{
		versormu[mu]=-1;	versormenosmu[mu]=1;
	}
	if(sentidonu==FRENTE){
		versornu[nu]=1;	versormenosnu[nu]=-1;
	}
	else{
		versornu[nu]=-1;	versormenosnu[nu]=1;
	}

	SomaVetoresPosicao(posicao,versormu,posicaomaismu);
	SomaVetoresPosicao(posicao,versornu,posicaomaisnu);
	SomaVetoresPosicao(posicaomaismu,versornu,posicaomaismumaisnu);
	SomaVetoresPosicao(posicaomaismu,versormenosnu,posicaomaismumenosnu);
	SomaVetoresPosicao(posicao,versormenosnu,posicaomenosnu);

}

void EloVizinhoSU3(int posicao[d],int direcao,int sentido,double u[3][3][2]){

// Calcula a matriz de elo vizinho a posição posicao na direção direcao e com sentido FRENTE ou TRAS

	double v[3][3][2];
	int posicaoaux[d];
	int versordir[d]={0,0,0,0};


	if(sentido==FRENTE) //	Elo no sentido positivo	é simplesmente o que está armazenado em U
	 	CopiarSU3(U[posicao[0]][posicao[1]][posicao[2]][posicao[3]][direcao],u);


	else{	//	Elo no sentido negativo é o conjugado de U da posição anterior
		versordir[direcao]=-1;	

		SomaVetoresPosicao(posicao,versordir,posicaoaux);
		ConjHermSU3(U[posicaoaux[0]][posicaoaux[1]][posicaoaux[2]][posicaoaux[3]][direcao],u);
	}
}