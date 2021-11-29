#ifndef FUNCOESINVERSAODIRAC_H
#define FUNCOESINVERSAODIRAC_H

void InicializarDirac();

void CalcularQtrevo(int posicao[d],int mu,int nu,double Q[3][3][2]);
void CalcularFtrevo(int posicao[d],double F[4][4][3][3][2]);

void CalcularcSWkappaSigmaF();

void CopiarVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);
void ProdutoEscalarVetorInversao(double num[2],double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2],double numf[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);
void SomaVetorInversao (double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);
void DiferencaVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double fmenosg[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);
void ProdutoVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double produto[2]);


void ProdutoDuploVetorInversao(double f1[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g1[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double produto1[2],double f2[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double g2[Nt][Nxyz][Nxyz][Nxyz][4][3][2], double produto2[2]);
void DiferencaVetorProdutoEscalarVetorInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double fmenosxg[Nt][Nxyz][Nxyz][Nxyz][4][3][2]);
void AcumularProdutoEscalarVetorInversao(double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2]);
void AcumularDuploProdutoEscalarVetorInversao(double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double y[2],double h[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2]);
void SomaProdutoVetorDiferencaProdutoInversao(double f[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double x[2],double g[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double y[2],double h[Nt][Nxyz][Nxyz][Nxyz][4][3][2],double i[Nt][Nxyz][Nxyz][Nxyz][4][3][2]);


void ProdutoDInversao(double f[Nt][Nxyz][Nxyz][Nxyz][d][3][2], double g[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);

void Inicializarx();
void Inicializarb(int inddirac, int indcor, double b[Nt][Nxyz][Nxyz][Nxyz][d][3][2]);

void BiCGStab(int inddirac, int indcor);

void ImprimirInverso(char nomearquivoinverso[ComprimentoMaxArquivo]);

void CarregarInverso(char nomearquivoinverso[ComprimentoMaxArquivo]);

void ImprimirOpDirac(char nomearquivoOpDirac[ComprimentoMaxArquivo]);

void ImprimirDelta(char nomearquivoinverso[ComprimentoMaxArquivo], double Delta[d][3][Nt][Nxyz][Nxyz][Nxyz][d][3][2]);


double ReTrDirac(double a[4][4][2]);
double ImTrProdutoDirac(double a[4][4][2],double b[4][4][2]);

void Calcularakslash(double amomento[d],double akslash[4][4][2]);
double akquadrado(double amomento[d]);

int cortemomentos(double ap[d]);

void CalcularFatoresdeForma(double amomento[d],double InversoMomento[4][4][2],char nomearquivoMZ[ComprimentoMaxArquivo]);
void TransformadaFourier(char nomearquivoMZ[ComprimentoMaxArquivo]);

#endif