//Parâmetros da rede

//Quantidade de pontos é 1 a mais que quantidade de elos

#define Nxyz 4 //Lado da rede na direção espacial
#define Nt 4 //Lado da rede na direção temporal
#define d 4	// Dimensão espaço temporal da rede

#define Volume Nxyz*Nxyz*Nxyz*Nt //	Número de pontos da rede

#define Nc 3.0	// Número de cores da teoria de gauge

#define FRIA 0	//Possíveis Partidas para Inicialização da Rede
#define QUENTE 1

#define TRAS -1
#define FRENTE 1

#define RE 0
#define IM 1

#define NAO 0
#define SIM 1

// Paramêtros gerais da Simulação

#define MAXconfigs 20 //	Número máximo de configurações descorrelacionadas a serem geradas


#define ComprimentoMaxArquivo 2000

#define kMAXbinomial 4