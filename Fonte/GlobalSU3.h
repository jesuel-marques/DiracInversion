#include "mtwister.h"
#include "../ParametrosSU3.h"
#include <time.h>

#ifndef GLOBALSU3_H
#define GLOBALSU3_H

extern double beta;

extern double UmComplexo[2];
extern double ZeroComplexo[2];

extern double identidadeSU2[4];
extern double nuloSU2[4];

extern double identidadeSU3[3][3][2];
extern double nuloSU3[3][3][2];

extern double R[3][3][2];
extern double S[3][3][2];
extern double T[3][3][2];

extern double U[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];
extern double Uaux[Nt][Nxyz][Nxyz][Nxyz][d][3][3][2];
extern double G[Nt][Nxyz][Nxyz][Nxyz][3][3][2];

//extern double tolerancia;

extern MTRand Aleatorio;

// extern double lambda[9][3][3][2];

#endif