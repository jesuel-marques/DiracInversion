#ifndef MEDICOESSU3_H
#define MEDICOESSU3_H

void SomaGramposSU3(int posicao[d],int mu,double A[3][3][2]);

void PlaquetasSU3(int posicao[d], int mu, int sentidomu, int nu, int sentidonu, double plaqueta[3][3][2]);
double MedirPlaquetasSU3();

void CalcularASU3(int posicao[d],int mu, double A[3][3][2]);
void DivergenciaASU3(int posicao[d], double divA[3][3][2]);


#endif