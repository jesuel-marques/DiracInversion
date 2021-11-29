#include <stdio.h>
#include <math.h>

#include "../ParametrosSU3.h"
#include "../ParametrosInversaoDirac.h"

#include "GlobalSU3.h"
#include "GlobalInversaoDirac.h"

#include "FuncoesMatematicas.h"
#include "SU3.h"
#include "RedeSU3.h"
#include "MedicoesSU3.h"


double Gama[d][4][4][2]; // Matrizes gama de Dirac
double Sigma[d][d][4][4][2];	//	Matrizes sigma de Dirac

double IdentidadeDirac[4][4][2]; // Identidade no espaço de Dirac

double IdentidadeDiracMenosGama[d][4][4][2];	//	Identidade menos gama de Dirac
double IdentidadeDiracMaisGama[d][4][4][2];	//	Identidade mais gama de Dirac


// Variáveis do algoritmo de inversão

double x[Nt][Nxyz][Nxyz][Nxyz][4][3][2];

double r[Nt][Nxyz][Nxyz][Nxyz][4][3][2]; 
double rtilde[Nt][Nxyz][Nxyz][Nxyz][4][3][2];

double v[Nt][Nxyz][Nxyz][Nxyz][4][3][2]; 
double p[Nt][Nxyz][Nxyz][Nxyz][4][3][2]; 

double aux1[Nt][Nxyz][Nxyz][Nxyz][4][3][2];	// b e t na notação do Leandro 
double aux2[Nt][Nxyz][Nxyz][Nxyz][4][3][2];	// s na notação do Leandro

void InicializarDirac(){

	//	Inicializa matrizes no espaço de Dirac (matrizes gamma e sigma)

	double comutador[4][4][2];	//	usado para calcular matrizes sigma
	double constantesigma[2]={0,-0.5};	//	usado para calcular matriz sigma

	int scan;

	//	Coloca os valores corretos nas matrizes gama de Dirac,
	//	assim como na identidade, identidade mais e menos gamas

	for(int alfa=0;alfa<4;alfa++)
		for(int beta=0;beta<4;beta++)
			for(int mu=0;mu<d;mu++){
				CopiarComplexo(ZeroComplexo,Gama[mu][alfa][beta]);
				for(int nu=0;nu<d;nu++)
					CopiarComplexo(ZeroComplexo,Sigma[mu][nu][alfa][beta]);
			}

	//	Seguindo convenções do Livro do Gattringer

	Gama[0][0][2][RE]=1.0;
	Gama[0][1][3][RE]=1.0;
	Gama[0][2][0][RE]=1.0;
	Gama[0][3][1][RE]=1.0;

	Gama[1][0][3][IM]=-1.0;
	Gama[1][1][2][IM]=-1.0;
	Gama[1][2][1][IM]=1.0;
	Gama[1][3][0][IM]=1.0;

	Gama[2][0][3][RE]=-1.0;
	Gama[2][1][2][RE]=1.0;
	Gama[2][2][1][RE]=1.0;
	Gama[2][3][0][RE]=-1.0;

	Gama[3][0][2][IM]=-1.0;
	Gama[3][1][3][IM]=1.0;
	Gama[3][2][0][IM]=1.0;
	Gama[3][3][1][IM]=-1.0;

	for(int alfa=0;alfa<4;alfa++)
		for(int beta=0;beta<4;beta++){
			if(alfa==beta)
				CopiarComplexo(UmComplexo,IdentidadeDirac[alfa][beta]);
		
			else
				CopiarComplexo(ZeroComplexo,IdentidadeDirac[alfa][beta]);


			for(int mu=0;mu<d;mu++){
				CopiarComplexo(IdentidadeDirac[alfa][beta],IdentidadeDiracMaisGama[mu][alfa][beta]);
				SomaComplexo(Gama[mu][alfa][beta],IdentidadeDiracMaisGama[mu][alfa][beta]);
				DiferencaComplexo(IdentidadeDirac[alfa][beta],Gama[mu][alfa][beta],IdentidadeDiracMenosGama[mu][alfa][beta]);
			}
		}

	//	Calculo das matrizes sigma a partir das matrizes gama
	//	Usando definição de Gattringer como sigma_munu=1/2i [gama_mu,gama_nu]
	//	A definição das matrizes gama de hep-lat/9605038 é diferente da
	//	que estou usando aqui, tirada de Gattringer.

	for(int mu=0;mu<d;mu++)
	for(int nu=0;nu<d;nu++){
		ComutadorComplexo4x4(Gama[mu],Gama[nu],comutador);
		ProdutocomEscalarComplexo4x4(constantesigma,comutador,Sigma[mu][nu]);
	}

}


void CalcularQtrevo(int posicao[d],int mu,int nu,double Q[3][3][2]){

	//	Usado para obter discretização de trevo do tensor de força

	double plaq1[3][3][2],plaq2[3][3][2],plaq3[3][3][2],plaq4[3][3][2];

	for(int a=0;a<3;a++)
	for(int b=0;b<3;b++)
		CopiarComplexo(ZeroComplexo,Q[a][b]);
	
	PlaquetasSU3(posicao,mu,FRENTE,nu,FRENTE,plaq1);
	PlaquetasSU3(posicao,nu,FRENTE,mu,TRAS,plaq2);
	PlaquetasSU3(posicao,mu,TRAS,nu,TRAS,plaq3);
	PlaquetasSU3(posicao,nu,TRAS,mu,FRENTE,plaq4);

	SomaSU3(plaq1,Q);
	SomaSU3(plaq2,Q);
	SomaSU3(plaq3,Q);
	SomaSU3(plaq4,Q);
}


void CalcularFtrevo(int posicao[d],double F[4][4][3][3][2]){

	//	Discretização em trevo do tensor de força.
	//	Calcula a matriz F que guarda o valor para todos os pontos e combinações mu, nu da rede.
	//	Calculado apenas para mu<nu por ser antisimétrico.
	//	Utilizando as convenções do livro de Gattringer.

	double constante[2]={0,-1.0/8.0};
	double Q1[3][3][2], Q2[3][3][2];
	double aux[3][3][2];
	
	for(int mu=0;mu<d;mu++)
	for(int nu=0;nu<d;nu++)
		if(mu<nu){

			CalcularQtrevo(posicao,mu,nu,Q1);
			CalcularQtrevo(posicao,nu,mu,Q2);

			DiferencaSU3(Q1,Q2,aux);
			MultiplicacaoEscalarSU3(aux,constante,F[mu][nu]);

		}


}

void CalcularcSWkappaSigmaF(){

	//	Termo de melhoria. Calculado imediatamente após carregar configuração e 
	//	mantido em memória para rápido acesso durante cálculo do algoritmo BiCGStab

	int a;
	double produto[2];
	int posicao[d];
	double F[4][4][3][3][2];	//	Tensor de força cromoeletromagnetico. A ser calculado a partir da discretização de trevo

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++){
		CalcularFtrevo(posicao,F);

		for(int alfa=0;alfa<4;alfa++)
		for(int alfal=0;alfal<4;alfal++)
			for(int a=0;a<3;a++)
			for(int al=0;al<3;al++){
				CopiarComplexo(ZeroComplexo,cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al]);
				for(int mu=0;mu<d;mu++)
				for(int nu=0;nu<d;nu++)
				if(mu<nu){
					ProdutoComplexoTres(cSWvezeskappa,Sigma[mu][nu][alfa][alfal],F[mu][nu][a][al],produto);
					SomaComplexo(produto,cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al]);
				}

			}
	}

}

// Aritmética dos vetores usados na inversão

void CopiarVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	// Copia cada entrada do vetor f para g

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				CopiarComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
}

void ProdutoEscalarVetorInversao(double num[2],double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double numf[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	// Realiza o produto num*f, isso é, multiplica cada entrada de f pelo numero complexo num, e coloca resultado em numf. 

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				ProdutoComplexo(num,f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],numf[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
}

void SomaVetorInversao (double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	//	Soma vetores f e g e acumula resultado em g.

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				SomaComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
}

void DiferencaVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double fmenosg[Nt][Nxyz][Nxyz][Nxyz][d][3][2]){
	
	//	Calcula a diferença entre cada entrada de f e g e coloca resultado em fmenosg

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				DiferencaComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],fmenosg[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
			
}

void ProdutoVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double produto[2]){
	
	//	Calcula o produto interno entre vetores. Isso será o resultado de fdagger.g, que é colocado na variável complexa produto

	int posicao[d];
	double termo[2];

	CopiarComplexo(ZeroComplexo,produto);

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
				ProdutoConjugadoComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],termo);
				SomaComplexo(termo,produto);
			}
}

void ProdutoDuploVetorInversao(double f1[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g1[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double produto1[2],double f2[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g2[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double produto2[2]){
	
	//	Calcula o produto interno entre vetores. Isso será o resultado de fdagger.g, que é colocado na variável complexa produto
	// Faz dois produtos independentes paralelamente

	int posicao[d];

	double termo[2];

	CopiarComplexo(ZeroComplexo,produto1);
	CopiarComplexo(ZeroComplexo,produto2);

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
				ProdutoConjugadoComplexo(f1[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g1[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],termo);
				SomaComplexo(termo,produto1);

				ProdutoConjugadoComplexo(f2[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g2[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],termo);
				SomaComplexo(termo,produto2);
			}
}

void DiferencaVetorProdutoEscalarVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double fmenosxg[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	//	fmenosxg=f-x*g
	//	calcula diferença entre vetor f e escalar complexo vezes vetor g
	
	int posicao[d];
	double prod[2];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
				ProdutoComplexo(x,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],prod);
				DiferencaComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],prod,fmenosxg[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
			}
}

void AcumularProdutoEscalarVetorInversao(double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	//	f=f+x*g
	//	Acumula resultado de multiplicação de escalar complexo x com vetor g em f
	
	int posicao[d];
	double prod[2];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){

				ProdutoComplexo(x,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],prod);
				SomaComplexo(prod,f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);

			}
}

void AcumularDuploProdutoEscalarVetorInversao(double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double y[2],double h[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	//	f=f+x*g+y*h
	//	Acumula resultado de multiplicação de escalar complexo x com vetor g
	//	e de escalar complexo y com vetor h em f

	int posicao[d];
	double prod1[2],prod2[2];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
				ProdutoComplexo(x,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],prod1);
				SomaComplexo(prod1,f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);

				ProdutoComplexo(y,h[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],prod2);
				SomaComplexo(prod2,f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
			}
}

void SomaProdutoVetorDiferencaProdutoInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double y[2],double h[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double i[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){
	
	// 	i=f+x*(g-y*h)

	int posicao[d];

	double yh[2];
	double termo[2];
	double xtermo[2];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
				CopiarComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],i[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);

				ProdutoComplexo(y,h[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],yh);
				DiferencaComplexo(g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],yh,termo);
				ProdutoComplexo(x,termo,xtermo);

				SomaComplexo(xtermo,i[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
			}
}

void ProdutoDInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){

	//	Calcula g=D.f, ou seja, a matriz de Dirac agindo no vetor f, e coloca o resultado em g.

	int posicao[d];
	int posicaol1[d];
	int posicaol2[d];

	int versormu[d];
	int versormenosmu[d];

	double Elo1[3][3][2];
	double produto1[2];

	double Elo2[3][3][2];
	double produto2[2];

	//	Para implementar condição antiperiódica

	double produtoaux[2];
	double MenosUmComplexo[2]={-1.0,0.0};

	for(int mu=0;mu<d;mu++){
		versormu[mu]=0;
		versormenosmu[mu]=0;
	}
							
	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){

				// Termo diagonal, simplesmente copia o conteúdo de f para g

				CopiarComplexo(f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);

				// Termo de hopping com vizinhos da "frente" (mu positivo, 1) e de "tras" (mu negativo, 2)

				for(int mu=0;mu<d;mu++){
					versormu[mu]=1;
					versormenosmu[mu]=-1;

					EloVizinhoSU3(posicao,mu,FRENTE,Elo1);
					EloVizinhoSU3(posicao,mu,TRAS,Elo2);

					SomaVetoresPosicao(posicao,versormu,posicaol1);
					SomaVetoresPosicao(posicao,versormenosmu,posicaol2);

					for(int alfal=0;alfal<4;alfal++)
						for(int al=0;al<3;al++){
							ProdutoComplexoQuatro(menoskappa,IdentidadeDiracMenosGama[mu][alfa][alfal],Elo1[a][al],f[posicaol1[0]][posicaol1[1]][posicaol1[2]][posicaol1[3]][alfal][al],produto1);
							ProdutoComplexoQuatro(menoskappa,IdentidadeDiracMaisGama[mu][alfa][alfal],Elo2[a][al],f[posicaol2[0]][posicaol2[1]][posicaol2[2]][posicaol2[3]][alfal][al],produto2);


							//	Condição de contorno antiperiódica na direção temporal
							if(posicao[0]!=(Nt-1)||mu!=0){
								SomaComplexo(produto1,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
							}
							else{
								ProdutoComplexo(produto1,MenosUmComplexo,produtoaux);
								SomaComplexo(produtoaux,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
							}

							if(posicao[0]!=0||mu!=0){
								SomaComplexo(produto2,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
							}
							else{
								ProdutoComplexo(produto2,MenosUmComplexo,produtoaux);
								SomaComplexo(produtoaux,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
							}
						}

					versormu[mu]=0;
					versormenosmu[mu]=0;
				}

				// Melhoria de trevo
		
				for(int alfal=0;alfal<4;alfal++)
					for(int al=0;al<3;al++){
						ProdutoComplexo(cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al],f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfal][al],produto1);
						SomaComplexo(produto1,g[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
				}
					

			}
}

void Inicializarx(){

	// Inicializa vetor x do algoritmo Bi-CGStab para inversão do operador de Dirac

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				for(int b=0;b<2;b++)

					//	Inicialização aleatória
					x[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][b]=genRand(&Aleatorio);
					
					//	Lüscher usa x_0=0 no openQCD. Mais rápido que inicalização aleatória
					// x[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][b]=0.0;
}	

void Inicializarb(int inddirac, int indcor, double b[Nt][Nxyz][Nxyz][Nxyz][4][3][2]){

	//	Inicializa b, que é um vetor coluna da matriz obtida pelas deltas de posição, índices de Dirac e cor

	int posicao[d];

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++)
				CopiarComplexo(ZeroComplexo,b[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a]);
			
	CopiarComplexo(UmComplexo,b[0][0][0][0][inddirac][indcor]);
}

void BiCGStab(int inddirac, int indcor){

	//	Algoritmo para inversão da matriz de Dirac (descrito em Gattringer pg. 140 e Leandro pg. 110)

	//	Algumas variáveis são globais por conta do tamanho necessário
	int cont=0;

	double auxa[2];
	double auxb[2];

	double alfa[2];
	double beta[2];
	double omega[2];

	double rho[2],rhoanterior[2];
	
	double norma;

	int primeiraiteracao=SIM;
	int invertido=NAO;

	Inicializarb(inddirac, indcor, aux1);
	
	// if(inddirac==0&&indcor==0)
		// Inicializarx();

	CopiarVetorInversao(aux1,x);/*x0=b*/

	ProdutoDInversao(x,aux2);/*aux2=A.x0*/
	DiferencaVetorInversao(aux1,aux2,r);/*rtilde=b-A.x0*/

	CopiarVetorInversao(r,rtilde);/*rtilde=r*/


	do{
		ProdutoVetorInversao(rtilde,r,rho);/*rho=rtildedag.r*/

		if(rho==0){
			printf("Metodo falhou.");
		}
		
		if(primeiraiteracao==SIM){
			CopiarVetorInversao(r,p);/*p=r*/
			primeiraiteracao=NAO;
		}

		else{
			ProdutoComplexo(alfa,rho,auxa);/*auxa=alfa*rho*/
			ProdutoComplexo(omega,rhoanterior,auxb);/*auxb=omega*rhoant*/
			DivisaoComplexo(auxa,auxb,beta);/*beta=alfa*rho/omega*rhoant*/

			SomaProdutoVetorDiferencaProdutoInversao(r,beta,p,omega,v,aux1);/*aux1=r+beta*(p-omega*v)*/
			CopiarVetorInversao(aux1,p);/*p=r+beta*(p-omega*v)*/
		}

		ProdutoDInversao(p,v);/*v=A.p*/
		ProdutoVetorInversao(rtilde,v,auxa);/*auxa=rtildedag.v*/
		DivisaoComplexo(rho,auxa,alfa);/*alfa=rho/(rtildedag.v)*/
		DiferencaVetorProdutoEscalarVetorInversao(r,alfa,v,aux2);/*s=r-alfa*v*/

		ProdutoVetorInversao(aux2,aux2,auxa);/*auxa=sdag.s*/
		norma=fabs(auxa[RE]);
		
		if(cont%5==0)
			printf("norma s: %.16lf\n",norma);

		if(norma<tolerancianormar){

			AcumularProdutoEscalarVetorInversao(alfa,p,x);/*x=x+alfa*p*/
			
			invertido=SIM;			
		}

		else{
			ProdutoDInversao(aux2,aux1);/*t=A.s*/

			ProdutoDuploVetorInversao(aux1,aux2,auxa,aux1,aux1,auxb);/*auxa=tdag.s; auxb=tdag.t*/
			DivisaoComplexo(auxa,auxb,omega);/*omega=tdag.t/tdag.s*/

			DiferencaVetorProdutoEscalarVetorInversao(aux2,omega,aux1,r);/*r=s-omega*t*/
			
			AcumularDuploProdutoEscalarVetorInversao(alfa,p,omega,aux2,x);/*x=x+alfa*p+omega*s*/

			CopiarComplexo(rho,rhoanterior);
		}

		cont++;

	}while(invertido==NAO);

	CopiarVetorInversao(x,Inverso[inddirac][indcor]);

}

void ImprimirInverso(char nomearquivoinverso[ComprimentoMaxArquivo]){

	//	Imprime em um arquivo as colunas da matriz inversa do operador de Dirac

	FILE *Arquivo;

	int posicao[d];

	Arquivo=fopen(nomearquivoinverso,"w+");

	//	Impressão em Arquivo binário

	// fwrite(Inverso, sizeof(Inverso),1, Arquivo);

	//	Impressão com índices e valores 

	// for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	// for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	// for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	// for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
	// 	for(int alfa=0;alfa<d;alfa++)
	// 	for(int beta=0;beta<d;beta++)
	// 		for(int a=0;a<3;a++)
	// 		for(int b=0;b<3;b++)
	// 			for(int c=0;c<2;c++)
	// 				fprintf(Arquivo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.16lf\n",posicao[0],posicao[1],posicao[2],posicao[3],alfa,beta,a,b,c,Inverso[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][c]);				

	//	Impressão em formato matricial (para comparação com Mathematica)

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++) //Linhas
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
			for(int a=0;a<3;a++){
					for(int beta=0;beta<4;beta++) //Colunas
					for(int b=0;b<3;b++){
				
					fprintf(Arquivo,"%.16lf+I(%.16lf)",Inverso[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][RE],Inverso[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][IM]);				
					fprintf(Arquivo,"\t");
					}
					fprintf(Arquivo,"\n");
				}
	fclose(Arquivo);
}

void CarregarInverso(char nomearquivoinverso[ComprimentoMaxArquivo]){

	//	Carrega, a partir de um arquivo, as colunas da matriz inversa do operador de Dirac

	FILE *Arquivo;

	int posicao[d];
	int alfa, beta;
	int a, b, c;
	double elementoleitura;

	Arquivo=fopen(nomearquivoinverso,"r");

	//	Lê de binário

	if(fread(Inverso,sizeof(Inverso),1,Arquivo)==1)
		printf("Arquivo de Inverso Carregado\n");

	//	Lê de arquivo com índices e valores
	// while(fscanf(Arquivo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",&posicao[0],&posicao[1],&posicao[2],&posicao[3],&alfa,&beta,&a,&b,&c,&elementoleitura)!=EOF)
	// 	Inverso[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][c]=elementoleitura;		
				
	fclose(Arquivo);
}

void ImprimirOpDirac(char nomearquivoOpDirac[ComprimentoMaxArquivo]){

	//	Imprime a Matriz de Dirac explicitamente (Comparação com Mathematica)

	FILE *Arquivo;

	int posicao[d];
	int posicao2[d];

	double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2];
	double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2];

	Arquivo=fopen(nomearquivoOpDirac,"w+");

	int cont1,cont2;


	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
	for(int alfa=0;alfa<4;alfa++)
	for(int a=0;a<3;a++){		

		f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][RE]=1.0;
		f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][IM]=0.0;
		ProdutoDInversao(f,g);
		f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][RE]=0.0;
		f[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][IM]=0.0;

		for(posicao2[0]=0;posicao2[0]<Nt;posicao2[0]++)
		for(posicao2[1]=0;posicao2[1]<Nxyz;posicao2[1]++)
		for(posicao2[2]=0;posicao2[2]<Nxyz;posicao2[2]++)
		for(posicao2[3]=0;posicao2[3]<Nxyz;posicao2[3]++)
		for(int beta=0;beta<4;beta++)
		for(int b=0;b<3;b++){
			fprintf(Arquivo,"%.20lf+I(%.20lf)",g[posicao2[0]][posicao2[1]][posicao2[2]][posicao2[3]][beta][b][RE],g[posicao2[0]][posicao2[1]][posicao2[2]][posicao2[3]][beta][b][IM]);
			fprintf(Arquivo,"\t");
		}
		fprintf(Arquivo,"\n");
	}		
	fclose(Arquivo);
}

void ImprimirDelta(char nomearquivoinverso[ComprimentoMaxArquivo], double Delta[d][3][Nt][Nxyz][Nxyz][Nxyz][d][3][2]){

	//	Imprime em um arquivo o produto entre a matriz de Dirac e a matriz inversa.
	//	Utilizado para testes da qualidade da inversão

	FILE *Arquivo;

	int posicao[d];

	double soma=0.0;

	Arquivo=fopen(nomearquivoinverso,"w+");

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
	for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
	for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
	for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
		for(int alfa=0;alfa<4;alfa++)
		for(int beta=0;beta<4;beta++)
			for(int a=0;a<3;a++)
			for(int b=0;b<3;b++)
				for(int c=0;c<2;c++){
					fprintf(Arquivo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.16lf\n",posicao[0],posicao[1],posicao[2],posicao[3],alfa,beta,a,b,c,Delta[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][c]);
					soma+=Delta[beta][b][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][c];
				}
	fclose(Arquivo);

	printf("%.16lf\n",soma);
}

double ReTrDirac(double a[4][4][2]){	

	//	Parte real do traço de uma matriz complexa 4x4
	
	return a[0][0][RE]+a[1][1][RE]+a[2][2][RE]+a[3][3][RE];
}

double ImTrProdutoDirac(double a[4][4][2],double b[4][4][2]){
	
	//	Parte imaginária do traço de um produto de matrizes complexas 4x4
	
	return a[0][0][RE]*b[0][0][IM] + a[0][0][IM]*b[0][0][RE] + 
		a[1][0][RE]*b[0][1][IM] + a[1][0][IM]*b[0][1][RE] + 
		a[2][0][RE]*b[0][2][IM] + a[2][0][IM]*b[0][2][RE] + 
		a[3][0][RE]*b[0][3][IM] + a[3][0][IM]*b[0][3][RE] + 
		a[0][1][RE]*b[1][0][IM] + a[0][1][IM]*b[1][0][RE] + 
		a[1][1][RE]*b[1][1][IM] + a[1][1][IM]*b[1][1][RE] + 
		a[2][1][RE]*b[1][2][IM] + a[2][1][IM]*b[1][2][RE] + 
		a[3][1][RE]*b[1][3][IM] + a[3][1][IM]*b[1][3][RE] + 
		a[0][2][RE]*b[2][0][IM] + a[0][2][IM]*b[2][0][RE] + 
		a[1][2][RE]*b[2][1][IM] + a[1][2][IM]*b[2][1][RE] + 
		a[2][2][RE]*b[2][2][IM] + a[2][2][IM]*b[2][2][RE] + 
		a[3][2][RE]*b[2][3][IM] + a[3][2][IM]*b[2][3][RE] + 
		a[0][3][RE]*b[3][0][IM] + a[0][3][IM]*b[3][0][RE] + 
		a[1][3][RE]*b[3][1][IM] + a[1][3][IM]*b[3][1][RE] + 
		a[2][3][RE]*b[3][2][IM] + a[2][3][IM]*b[3][2][RE] + 
		a[3][3][RE]*b[3][3][IM] + a[3][3][IM]*b[3][3][RE];
}



void Calcularakslash(double amomento[d],double akslash[4][4][2]){	
	
	// Calcula akslash=sum_mu ak_mu*gamma_mu a partir de ap
	
	double ak[d];

	for(int i=0;i<4;i++){
		ak[i]=sin(amomento[i]);
		for(int j=0;j<4;j++)
			CopiarComplexo(ZeroComplexo,akslash[i][j]);
	}

	akslash[0][2][RE]=ak[0];
	akslash[1][3][RE]=ak[0];
	akslash[2][0][RE]=ak[0];
	akslash[3][1][RE]=ak[0];

	akslash[0][3][IM]=-ak[1];
	akslash[1][2][IM]=-ak[1];
	akslash[2][1][IM]=ak[1];
	akslash[3][0][IM]=ak[1];

	akslash[0][3][RE]=-ak[2];
	akslash[1][2][RE]=ak[2];
	akslash[2][1][RE]=ak[2];
	akslash[3][0][RE]=-ak[2];

	akslash[0][2][IM]=-ak[3];
	akslash[1][3][IM]=ak[3];
	akslash[2][0][IM]=ak[3];
	akslash[3][1][IM]=-ak[3];
}

double akquadrado(double amomento[d]){	
	
	//	Calcula (ak)²=sum_mu (sin(ap_mu)^2)
	
	double ak2=0.0;
	for(int mu=0;mu<d;mu++){
		ak2+=pow2(sin(amomento[mu]));
	}
	return ak2;
}

int cortemomentos(double ap[d]){

	//	Condição de corte cilindrico que restringe os momentos àqueles próximos à diagonal
	//	Condição definida em hep-lat/0007028

	double n[d]={0.5,0.5,0.5,0.5};
	double angulocorte=M_PI/8.0;
	
	double ap2=0.0;
	double norma;
	double apn=0.0;

	for(int i=0;i<d;i++){
		apn+=ap[i]*n[i];
		ap2+=pow2(ap[i]);
	}
	norma=sqrt(ap2);

	if(norma*sin(acos(apn/norma))<=angulocorte)
		return 1;
	else
		return 0;
}

void CalcularFatoresdeForma(double amomento[d],double PropagadorQuarkMomento[4][4][2],char nomearquivoMZ[ComprimentoMaxArquivo]){

	//	Calcula os fatores de forma Z e M como definidos em hep-lat/0007028

	//	Arquivo no qual será impresso Z e M
	FILE *MeZ;
	//	Variáveis de momento usadas no cálculo. Definidas em hep-lat/0007028
	double akslash[4][4][2];
	double ak2;

	double calA,calB,A,B,Z,M;	//	cal se refere a caligrafico, como denotado no artigo
	double calAI,calBI,AI,BI,ZI,MI;


	double PropagadorQuarkMomentoMelhorado[4][4][2];	//	Propagador melhorado I do artigo. (1+am)S_0-1/2
	double MenosMeioComplexo[2]={-0.5,0.0};		//	Usado no cálculo do propagador melhorado


	for(int a=0;a<4;a++)
	for(int b=0;b<4;b++){
		for(int c=0;c<2;c++){
			PropagadorQuarkMomentoMelhorado[a][b][c]=(1.0+am)*PropagadorQuarkMomento[a][b][c];
		}

		if(a==b){
			SomaComplexo(MenosMeioComplexo,PropagadorQuarkMomentoMelhorado[a][b]);
		}
	}
	
	//	Determina valores de ak e outros relacionados, a partir de ap 
	
	Calcularakslash(amomento,akslash);
	ak2=akquadrado(amomento);

	//	Fatores de forma do propagador ingênuo

	calA=-ImTrProdutoDirac(akslash,PropagadorQuarkMomento)/(4.0*ak2);
	calB=ReTrDirac(PropagadorQuarkMomento)/(4.0);

	A=calA/(ak2*pow2(calA)+pow2(calB));
	B=calB/(ak2*pow2(calA)+pow2(calB));


	Z=1.0/A;
	M=B/A;

	//	Fatores de forma do propagador melhorado

	calAI=-ImTrProdutoDirac(akslash,PropagadorQuarkMomentoMelhorado)/(4.0*ak2);
	calBI=ReTrDirac(PropagadorQuarkMomentoMelhorado)/(4.0);

	AI=calAI/(ak2*pow2(calAI)+pow2(calBI));
	BI=calBI/(ak2*pow2(calAI)+pow2(calBI));


	ZI=1.0/AI;
	MI=BI/AI;
	
	//	Imprime valores do momento e dos fatores de forma
	MeZ=fopen(nomearquivoMZ,"a+");

	fprintf(MeZ,"%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",amomento[0],amomento[1],amomento[2],amomento[3],Z,M,ZI,MI);

	fclose(MeZ);

}

void TransformadaFourier(char nomearquivoMZ[ComprimentoMaxArquivo]){

	//	Realiza a transformada de Fourier do propagador no espaço x para o espaço p

	int posicao[d];
	double amomento[d];

	//	Para cada momento, o propagador será uma matriz 4x4 complexa

	double InversoMomento[4][4][2];
	double PropagadorQuarkMomento[4][4][2];	
	//	O propagador de quark é dado pela inversa com uma normalização.
	//	A normalização é necessária porque a matriz a ser invertida é dividida por am0+4.
	//	Isso é feito com o intuito de invertemos uma matriz tipo 1-kappa*hopping-kappa csw SW
	//	Além disso, como o traço no espaço de cor já está sendo realizado nessa função,
	//	também multiplico o resultado por 

	//	Variaveis da normalização

	double am0=1.0/(-2.0*menoskappa[RE])-4.0;

	double normalizacao[2]={(am0+4.0)*3.0,0.0};

	//	Variáveis auxiliares para cálculo da transformada
	double termo[2]={0.0,0.0};
	double expipx[2]={0.0,0.0};
	
	for(int nt=1;nt<=Nt;nt++){
    amomento[0]=2.0*((double)M_PI/Nt)*((double)nt-0.5-(double)Nt/2.0);	//	São calculados todos os momentos possíveis
    for(int nx=1;nx<=Nxyz;nx++){										//	permitidos pelas condições de contorno.
    amomento[1]=2.0*((double)M_PI/Nxyz)*((double)nx-(double)Nxyz/2.0);	//	As fórmulas são dadas em hep-lat/0007028
	for(int ny=1;ny<=Nxyz;ny++){
	amomento[2]=2.0*((double)M_PI/Nxyz)*((double)ny-(double)Nxyz/2.0);
    for(int nz=1;nz<=Nxyz;nz++){
    amomento[3]=2.0*((double)M_PI/Nxyz)*((double)nz-(double)Nxyz/2.0);

		// if(cortemomentos(amomento)){	//	Restrição aos momentos que obedecem a condição do corte cilíndrico
		
			for(int alfa=0;alfa<4;alfa++)
			for(int beta=0;beta<4;beta++)
				CopiarComplexo(ZeroComplexo,InversoMomento[alfa][beta]);
			
			//	Para cada momento faz a transformada
					
			for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)	
			for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
			for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
			for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)	

				for(int alfa=0;alfa<4;alfa++)
				for(int beta=0;beta<4;beta++)
				for(int a=0;a<3;a++){	//	Já tira o traço em cor nessa etapa.

					Calcularexp(amomento,posicao,expipx);	//	Calcula a exponencial

					
					ProdutoComplexo(Inverso[beta][a][posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a],expipx,termo);
					SomaComplexo(termo,InversoMomento[alfa][beta]); 

				}

				//	Normalização

				for(int alfa=0;alfa<4;alfa++)
				for(int beta=0;beta<4;beta++)
					DivisaoComplexo(InversoMomento[alfa][beta],normalizacao,PropagadorQuarkMomento[alfa][beta]);

				//	Calculo de fatores de forma Z e M e impressão deles num arquivo

				CalcularFatoresdeForma(amomento,PropagadorQuarkMomento,nomearquivoMZ);
				
		// }					
	}
    }
	}
    }
	
}