#include <iostream>
#include <armadillo>
#include <functions.h>

using namespace std;
using namespace arma;

void ForwardEuler();
void BackwardEuler();
void metropolis();

int main()
{
    int Nx, Nt, Ny; double dt, dx, x0, y0, xl, yl, t0; vec x, y, t;

    initialize(&Nx,&Ny,&Nt,&dt,&dx,&x,&y,&t);
    conditions(Nx,Ny,&x0, &y0, &xl, &yl, &t0, &x, &y, &t);



    // First attempt at Markov Chain
    int N0; double l0, r; vec X;

    N0 = 1000; X = zeros<vec>(N0); l0 = sqrt(2*dt);

    for (int t=0;t<=Nt;t++){
        for (int i=0;i<N0;i++){
            r = 0.5; // Insert random number generator

            if (X(i)==0){
                if (r<0.5) {X(i) += l0;}
                else       {}// Remove array element
                //Add new element with value 0
            }
            else if (X(i) > 1){
                if (r<0.5) {}// Remove array element
                else       {X(i) -= l0;}
            }
            else{
                if (r<0.5) {X(i) += l0;}
                else       {X(i) -= l0;}
            }
        }
    }


    return 0;
}


void ForwardEuler(){

}

