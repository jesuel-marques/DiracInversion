#ifndef SU3_H
#define SU3_H

void CopiarSU3(double u[3][3][2], double v[3][3][2]);

void SomaSU3(double u[3][3][2],double resultado[3][3][2]);
void DiferencaSU3(double u[3][3][2], double v[3][3][2], double umenosv[3][3][2]);

void TrSU3(double u[3][3][2],double tr[2]);
double ReTr(double u[3][3][2]);

void DeterminanteSU3(double u[3][3][2], double det[2]);

void ConjHermSU3(double X[3][3][2],double Xch[3][3][2]);

void ProdutoSU3(double Y[3][3][2], double Z[3][3][2], double YZ[3][3][2]);
void ProdutoSU3Tres(double X[3][3][2], double Y[3][3][2], double Z[3][3][2], double XYZ[3][3][2]);
void ProdutoSU3Quatro(double W[3][3][2], double X[3][3][2], double Y[3][3][2], double Z[3][3][2], double WXYZ[3][3][2]);

void PotenciaSU3(double u[3][3][2],double unapotencias[kMAXbinomial+1][3][3][2]);

void MultiplicacaoEscalarSU3(double u[3][3][2],double num[2], double numvezesu[3][3][2]);

void ProjecaoSU3(double u[3][3][2], double uSU3[3][3][2]);

void DecomporAlgebraSU3(double g[3][3][2],double ga[9]);

// void InicializarLambda();

void ImprimirMatrizSU3(double M[3][3][2],char *nome );

#endif