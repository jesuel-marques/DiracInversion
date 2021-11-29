#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "../ParametrosSU3.h"
#include "GlobalSU3.h"
#include "GlobalInversaoDirac.h"

#include "FuncoesMatematicas.h"

#include "RedeSU3.h"
#include "SU3.h"


void SomaGramposSU3(int posicao[d],int mu,double A[3][3][2]){

	//	Calcula soma dos grampos para calculo da ação no programa que gera as configurações.
	//	Dado um elo por posição e índice de espaço-tempo, devolve a soma das grampos em A 
	
	int posicaomaismu[d];
	int posicaomaismumaisnu[d];
	int posicaomaisnu[d];
	int posicaomaismumenosnu[d];
	int posicaomenosnu[d];


	double u1a[3][3][2], u1b[3][3][2], u1c[3][3][2], u1[3][3][2];
	double u2a[3][3][2], u2b[3][3][2],u2c[3][3][2], u2[3][3][2];

	CopiarSU3(nuloSU3,A);

	for(int nu=0;nu<d;nu++){
		if(mu!=nu){
			PosicoesVizinhas(posicao,mu,FRENTE,nu,FRENTE,posicaomaismu,posicaomaismumaisnu,posicaomaisnu,posicaomaismumenosnu,posicaomenosnu);

			EloVizinhoSU3(posicaomaismu,nu,FRENTE,u1a);
			EloVizinhoSU3(posicaomaismumaisnu,mu,TRAS,u1b);
			EloVizinhoSU3(posicaomaisnu,nu,TRAS,u1c);
			ProdutoSU3Tres(u1a,u1b,u1c,u1);

			SomaSU3(u1,A);

			EloVizinhoSU3(posicaomaismu,nu,TRAS,u2a);
			EloVizinhoSU3(posicaomaismumenosnu,mu,TRAS,u2b);
			EloVizinhoSU3(posicaomenosnu,nu,FRENTE,u2c);
			ProdutoSU3Tres(u2a,u2b,u2c,u2);
			
			SomaSU3(u2,A);
		}
	}
}

void PlaquetasSU3(int posicao[d], int mu, int sentidomu, int nu, int sentidonu, double plaqueta[3][3][2]){

	//	Dada uma posição e dois índices mu e nu (com seus sentidos podendo ser FRENTE e TRAS)
	//	devolve a plaqueta correspondente. Especificar o sentido é necessário para a
	//	discretização de trevo do tensor de força.
	
	int posicaomaismu[d];
	int posicaomaismumaisnu[d];
	int posicaomaisnu[d];
	int posicaomaismumenosnu[d];
	int posicaomenosnu[d];

	int a;

	int sentidomuinv,sentidonuinv;
	double ua[3][3][2], ub[3][3][2], uc[3][3][2], ud[3][3][2];

	if(sentidomu==FRENTE)
		sentidomuinv=TRAS;
	else
		sentidomuinv=FRENTE;

	if(sentidonu==FRENTE)
		sentidonuinv=TRAS;
	else
		sentidonuinv=FRENTE;
	

	PosicoesVizinhas(posicao,mu, sentidomu, nu,sentidonu, posicaomaismu,posicaomaismumaisnu,posicaomaisnu,posicaomaismumenosnu,posicaomenosnu);
	
	EloVizinhoSU3(posicao,mu,sentidomu,ua);
	EloVizinhoSU3(posicaomaismu,nu,sentidonu,ub);
	EloVizinhoSU3(posicaomaismumaisnu,mu,sentidomuinv,uc);
	EloVizinhoSU3(posicaomaisnu,nu,sentidonuinv,ud);

	ProdutoSU3Quatro(ua,ub,uc,ud,plaqueta);
}

double MedirPlaquetasSU3(){

//	Função que mede a soma sobre toda a rede da parte real do traço das plaquetas dividido por Nc=3.

	int posicao[d];

	double plaqueta[3][3][2];
	double traco[2];

	double somaretrplaquetas=0;

	int cont=0;

	for(posicao[0]=0;posicao[0]<Nt;posicao[0]++)
		for(posicao[1]=0;posicao[1]<Nxyz;posicao[1]++)
			for(posicao[2]=0;posicao[2]<Nxyz;posicao[2]++)
				for(posicao[3]=0;posicao[3]<Nxyz;posicao[3]++)
					for(int mu=0;mu<d;mu++)
						for(int nu=0;nu<d;nu++)
							if(mu<nu){

								PlaquetasSU3(posicao,mu,FRENTE,nu,FRENTE,plaqueta); // Função retorna plaqueta com posição dada e direções mu e nu 

								TrSU3(plaqueta,traco);
																
								somaretrplaquetas+=(double)traco[0]/Nc;

							}
							
	return somaretrplaquetas;

}


void CalcularASU3(int posicao[d],int mu, double A[3][3][2]){

//	Função que calcula o campo de gauge a partir dos elos em Uaux, que é a matriz transformada de gauge do programa de fixação.

	double Udag[3][3][2];
	double UmenosUdag[3][3][2];
	double escalar1[2]={0.0,-0.5};
	double escalar2[2]={-(double)1.0/3.0,0.0};

	double traco[2];
	double fator[2];
	double partetraco[3][3][2];

	ConjHermSU3(Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],Udag);
	DiferencaSU3(Uaux[posicao[0]][posicao[1]][posicao[2]][posicao[3]][mu],Udag,UmenosUdag);
	MultiplicacaoEscalarSU3(UmenosUdag,escalar1,A);

	//Tirar parte do traço

	TrSU3(A,traco);

	ProdutoComplexo(traco,escalar2,fator);
	MultiplicacaoEscalarSU3(identidadeSU3,fator,partetraco);// isso pode ser melhorado, toda vez calculando identidade vezes 1/Nc
	SomaSU3(partetraco,A);


}

void DivergenciaASU3(int posicao[d], double divA[3][3][2]){	
	
	//	Calcula a divergência na rede do campo de gauge.
	
	int versormenosmu[d]={0,0,0,0};
	int posicaomenosmu[d];

	double A1[3][3][2];	//	Calcula através da diferença entre As de elos separados por um elo. 
	double A2[3][3][2];

	double termodivA[3][3][2];

	CopiarSU3(nuloSU3,divA);


	for(int mu=0;mu<d;mu++){

		versormenosmu[mu]=-1;
		SomaVetoresPosicao(posicao,versormenosmu,posicaomenosmu);
		
		CalcularASU3(posicao,mu,A1);		//melhorar esse cálculo... estou calculando A muitas vezes, melhor fazer varredura calculando A da rede inteira e depois só ler
		CalcularASU3(posicaomenosmu,mu,A2);

		DiferencaSU3(A1,A2,termodivA);

		SomaSU3(termodivA,divA);	//	Soma os termos correspondendo as diferentes direções mu.

		versormenosmu[mu]=0;
	}

}