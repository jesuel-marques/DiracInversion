#ifndef REDESU3_H
#define REDESU3_H

void InicializarUSU3(int Partida);

void TransformacaoCalibreRandomicaSU3();

void CopiarConfigSU3(double ConfigU[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2],double ConfigUcopiado[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2]);

void CarregarUSU3(char NomeArquivo[ComprimentoMaxArquivo]);

void CarregarUSU3Matheus(char NomeArquivo[ComprimentoMaxArquivo]);

void ImprimirUSU3(char NomeArquivo[ComprimentoMaxArquivo]);

void CopiarPosicao(int v[d],int u[d]);

void SomaVetoresPosicao(int v[d], int u[d], int vmaisu[d]);

void PosicoesVizinhas(int posicao[d], int mu, int sentidomu, int nu, int sentidonu, int posicaomaismu[d], int posicaomaismumaisnu[d],int posicaomaisnu[d], int posicaomaismumenosnu[d], int posicaomenosnu[d]);

void EloVizinhoSU3(int posicao[d],int direcao,int sentido,double u[3][3][2]);

#endif