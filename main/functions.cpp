#include "functions.h"

void initialize(int *Nx, int *Nt, double *T){
    cout << "Nx: "; cin >> *Nx;
    cout << "Nt: "; cin >> *Nt;
    cout << "T : "; cin >> *T ;
}


void ForwardEuler(mat *U, int Nx, int Nt, double h, double dt){
    mat u, unew; double a;
    u = *U;
    unew = u;
    a = dt/h/h;

    for (int t = 0; t<= Nt; t++){
        for (int i = 1; i < Nx-1; i++){
            for (int j = 1; j < Nx-1; j++){
                unew(i,j) = u(i,j) + a*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
            }
        }
    u = unew;
    }
    *U = u;
}


void BackwardEuler(mat *U, int Nx, double h, double dt){
    mat unew, u; double a,a2,difference,tolerance; int n, max_iter;
    u = *U;
    a = dt/h/h; a2 = 1 / (1 + 4*a);
    n = 0; difference = 1; tolerance = 1e-6; max_iter = 1e6;

    while (n < max_iter && difference > tolerance){
        unew = u;
        difference = 0;
        for (int i=1; i<Nx-1; i++){
            for (int j=1; j<Nx-1; j++){
                unew(i,j) = a2*(u(i,j) + a*(unew(i+1,j) + unew(i-1,j) + unew(i,j+1) + unew(i,j-1)));
                difference += abs(unew(i,j)-u(i,j));
            }
        }
        difference = difference / Nx/Nx;
        n ++;
        u = unew;
    }
    *U = u;
}
