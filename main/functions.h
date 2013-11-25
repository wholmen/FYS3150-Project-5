#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#endif // FUNCTIONS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void initialize(int *Nx, int *Ny, int *Nt, double *dt, double *dx, vec *x, vec *y, vec *t);
void conditions(int Nx, int Ny, double *x0, double *y0, double *xl, double *yl, double *t0, vec *x, vec *y, vec *t);
