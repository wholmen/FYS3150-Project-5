#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#endif // FUNCTIONS_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void initialize(int *Nx, int *Nt, double *T);
void ForwardEuler(mat *U, int Nx, int Nt, double h, double dt);
void BackwardEuler(mat *U, int Nx, double h, double dt);
