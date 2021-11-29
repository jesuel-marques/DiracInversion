#ifndef GLOBALINVERSAODIRAC_H
#define GLOBALINVERSAODIRAC_H

extern double menoskappa[2];	//	Parâmetro associado a massa dos fermions de Wilson. Gattringer Eq. 5.55 
extern double cSWvezeskappa[2];
extern double am;


extern double Inverso[4][3][Nt][Nxyz][Nxyz][Nxyz][4][3][2];	// Matriz que conterá as colunas relevantes do inverso da matriz de Dirac
extern double cSWkappaSigmaF[Nt][Nxyz][Nxyz][Nxyz][4][3][4][3][2]; //soma_munu sigma_munu F_munu 
#endif