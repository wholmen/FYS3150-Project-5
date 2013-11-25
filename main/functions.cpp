#include "functions.h"

void initialize(int *Nx, int *Ny, int *Nt, double *dt, double *dx, vec *x, vec *y, vec *t){
    cout << "Nx: "; cin >> *Nx;
    cout << "Ny: "; cin >> *Ny;
    cout << "Nt: "; cin >> *Nt;
    cout << "dt: "; cin >> *dt;
    cout << "dx: "; cin >> *dx;
    *x = zeros<vec>(*Nx); *y = zeros<vec>(*Ny); *t = zeros<vec>(*Nt);
}

void conditions(int Nx, int Ny, double *x0, double *y0, double *xl, double *yl, double *t0, vec *x, vec *y, vec *t){
    int test; vec X = *x; vec Y = *y; vec T = *t;
    cout << "use default conditions? Yes: 0. No: 1"; cin >> test;
    if (test == 1){
        cout << "x0: "; cin >> *x0;
        cout << "y0: "; cin >> *y0;
        cout << "xl: "; cin >> *xl;
        cout << "yl: "; cin >> *yl;
        cout << "t0: "; cin >> *t0;
    }
    else{
        *x0 = 1; *y0 = 1; *xl = 0; *yl = 0; *t0 = 0;
    }
    X(0) = *x0; X(Nx-1) = *xl;
    Y(0) = *y0; Y(Ny-1) = *yl;
    T(0) = *t0;
    *x = X, *y = Y; *t = T;
}
