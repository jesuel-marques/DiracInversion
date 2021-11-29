//gcc -o InversaoDiracInteragente Fonte/RedeSU3.c Fonte/FuncoesMatematicas.c Fonte/mtwister.c Fonte/SU3.c Fonte/FuncoesInversaoDirac.c InversaoDiracInteragente.c Fonte/ranlxd.c Fonte/ranlux_common.c Fonte/MedicoesSU3.c -lm -O3

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ParametrosSU3.h"
#include "ParametrosInversaoDirac.h"

#include "Fonte/GlobalSU3.h"
#include "Fonte/GlobalInversaoDirac.h"

#include "Fonte/mtwister.h"
#include "Fonte/ranlux.h"

#include "Fonte/FuncoesInversaoDirac.h"

#include "Fonte/RedeSU3.h"
#include "Fonte/MedicoesSU3.h"

double U[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];	//	Matriz de elos, utilizados para calcular matriz de Dirac
double Uaux[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];	//	Matrix com o links atualizados pela transformação de calibre. 
double G[Nt][Nxyz][Nxyz][Nxyz][3][3][2];	//	Transformação de calibre sobre U(x) para levar ao calibre de Landau

double cSWkappaSigmaF[Nt][Nxyz][Nxyz][Nxyz][4][3][4][3][2]; //soma_munu sigma_munu F_munu 

double beta=5.9;	//	Inverso do acoplamento

double Inverso[4][3][Nt][Nxyz][Nxyz][Nxyz][4][3][2]; // Matriz que conterá as colunas relevantes do inverso da matriz de Dirac

MTRand Aleatorio;	//	Variável aleatória do gerador MersenneTwister

//	Parâmetro associado a massa dos fermions de Wilson. Gattringer Eq. 5.55 

// double menoskappa[2]={-0.12314361007807304879,0.0};	//	0.123... é usado para se obter massa nua no caso livre de am=0.0603. 
double menoskappa[2]={-0.1369017589,0.0};				//  0.1369... usado para obter am=0.0603
//  double menoskappa[2]={-0.137,0.0};					//	0.137... valor dado por Jon-Ivar no artigo hep-lat/0007028. Resulta em am=0.057681

// Coeficiente de Sheikholeslami-Wohlert vezes kappa. É o que aparece na ação
double cSWvezeskappa[2];

//	Massa em unidades da rede. É calculada a partir de 1/(2 kappa)-1/(2 kappa_c). 
//	Para teoria livre kappa_c=0.125. Para teoria interagente em hep-lat/0007028, kappa_c=0.1392

double am;

double Delta[4][3][Nt][Nxyz][Nxyz][Nxyz][4][3][2];

int main(){

	char nomearquivoConfig[ComprimentoMaxArquivo-200];	//	Nome do arquivo de configuração de elos a ser carregada
	char nomearquivoConfigCompleto[ComprimentoMaxArquivo];	//	Nome do arquivo de configuração de elos com pasta

	char nomearquivoOpDirac[ComprimentoMaxArquivo];	//	Nome do arquivo a conter o operador de Dirac 
	char nomearquivoinverso[ComprimentoMaxArquivo];	//	Nome do arquivo a conter a inversa da matriz de Dirac

	char nomearquivodelta[ComprimentoMaxArquivo];

	char nomearquivoMZ[ComprimentoMaxArquivo];


	//	Inicialização de geradores de números aleatórios	
	Aleatorio=seedRand(time(NULL));	//	Inicialização do gerador Aleatório MersenneTwister
	rlxd_init(1,97021123);


	//	Massa em unidades da rede. É calculada a partir de 1/(2 kappa)-1/(2 kappa_c). 
	//	Para teoria livre kappa_c=0.125. Para teoria interagente em hep-lat/0007028, kappa_c=0.1392


	am=1.0/(-2.0*menoskappa[RE])-1.0/(2.0*0.1392);
	// am=1.0/(-2.0*menoskappa[RE])-1.0/(2.0*0.125);

	//	Termo de Sheikholeslami-Wohlert vezes kappa é o que aparece na ação melhorada.
	//	Valor de c_SW=1.479 é dado em hep-lat/0007028

	cSWvezeskappa[RE]=1.479*(-menoskappa[RE]);
	cSWvezeskappa[IM]=0;

	//	Inicializa as matrizes do espaço de Dirac.

	InicializarDirac();	

	printf("kappa: %lf, a*m: %lf\n",-menoskappa[0],am);

	//	Carrega configuração U trivial igual a identidade.

	InicializarUSU3(FRIA);

	for(int config=1;config<=2;config++){

		//	Gera nomes das configurações a serem carregadas
		
		sprintf(nomearquivoConfig,"Config_%d_beta_%.3lf_Nxyz_%d_Nt_%d.txt",config,beta,Nxyz,Nt);
		sprintf(nomearquivoConfigCompleto,"./Configs/%s",nomearquivoConfig);

		//	Carrega as diferentes configurações de elos para gerar matriz de Dirac a ser invertida

		CarregarUSU3(nomearquivoConfigCompleto);
	
		//	Tendo lido a configuração, pode-se calcular o tensor de força

		CalcularcSWkappaSigmaF();

		sprintf(nomearquivoOpDirac,"./DiracOperators/OpDiracInteracting_config_%d_beta_%.3lf_Nxyz_%d_Nt_%d_kappa_%.6lf.txt",config,beta,Nxyz,Nt,-menoskappa[RE]);

		ImprimirOpDirac(nomearquivoOpDirac);	/*Imprime Operador de Dirac para testes (Comparação com Mathematica)*/

		//	Produz as doze (4(Dirac)*3(Cor)) colunas da matriz inversa de Dirac

		for(int inddirac=0;inddirac<4;inddirac++)
			for(int indcor=0;indcor<3;indcor++){
				printf("Dirac: %d, Cor: %d\n",inddirac,indcor);
				BiCGStab(inddirac,indcor);	// Para cada combinação de índice spinorial e de cor, produz coluna da matriz inversa do Operador de Dirac
			}

		//	Gera nome e imprime resultado da inversão da matriz de Dirac (Comparação com Mathematica)

		sprintf(nomearquivoinverso,"./DiracInverses/InverseInteracting_config_%d_beta_%.3lf_Nxyz_%d_Nt_%d_kappa_%.6lf.txt",config,beta,Nxyz,Nt,-menoskappa[RE]);

		ImprimirInverso(nomearquivoinverso); /* Imprime Inverso do Operador de Dirac para testes (Comparação com Mathematica) */

		for(int inddirac=0;inddirac<4;inddirac++)
			for(int indcor=0;indcor<3;indcor++){
				ProdutoDInversao(Inverso[inddirac][indcor],Delta[inddirac][indcor]);
			}

		sprintf(nomearquivodelta,"./Deltas/DeltaInteracting_config_%d_beta_%.3lf_Nxyz_%d_Nt_%d_kappa_%.6lf.txt",config,beta,Nxyz,Nt,-menoskappa[RE]);
		ImprimirDelta(nomearquivodelta,Delta);


		sprintf(nomearquivoMZ,"./MZ/MZInteracting_config_%d_Nxyz_%d_Nt_%d_kappa_%.6lf.txt",config,Nxyz,Nt,-menoskappa[RE]);

		TransformadaFourier(nomearquivoMZ);

	}

	return 0;
}