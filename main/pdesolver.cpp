#include "pdesolver.h"

PDESolver::PDESolver(int nt, double Dt)
{
    Nt = nt; dt = Dt;
}


void PDESolver::RandomWalk(vec *x)
{
    double l0, r; vec X; long idum; int N,n;
    X = *x; l0 = sqrt(2*dt); idum = -1;

    for (int t=0;t<=Nt;t++){
        N = X.n_rows; n = 0;            // n will determine how many new particles to be added in x = 0

        for (int i=0;i<N;i++){
            r = ran0(&idum);            // Getting a random number

            if (X(i) < 1e-10) {n += 1;} // This particle will move away from x0. Need to add one more.

            if (r>0.5) {X(i) += l0;}    // Move to the right
            else       {X(i) -= l0;}    // Move to the left

            if (abs(X(i)) < 1e-10) {n -= 1;} // This particle has arrived at x0. Need to add one less
        }
        int i = 0;
        while (i<N){
            if      (X(i) > 1){X.shed_row(i); N-=1;} // Remove particle and reduce total particle #
            else if (X(i) < 0){X.shed_row(i); N-=1;} // Remove particle and reduce total particle #
            else              {i++;}
        }
        for (int i=0;i<n;i++){
            X.insert_rows(0,1); X(0) = 0; // Maintaining # of particles at x=0
        }
    }
    *x = sort(X);
}
