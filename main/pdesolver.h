#ifndef PDESOLVER_H
#define PDESOLVER_H
#include <iostream>
#include <armadillo>
#include "lib.h"

using namespace std;
using namespace arma;

class PDESolver
{
    int Nt; double dt;

public:
    PDESolver(int nt, double Dt);
    void RandomWalk(vec *x);
};

#endif // PDESOLVER_H
