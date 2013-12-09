#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <iostream>
#include <armadillo>
#include "lib.h"

using namespace std;
using namespace arma;

class PDESolver
{
    int Nt; double dt; int Nx; double dx; double pi = 3.14159265359;

public:
    PDESolver(int nt, int nx, double Dt, double Dx);
    void RandomWalk(vec *x);
    void GaussRandomWalk(vec *x);
    void Analytical(vec *U);
    void CrankNicolson(vec *v);
    vec tridiagonal(vec a1, vec a2, vec a3, vec b);
};

#endif // PDESOLVER_H
