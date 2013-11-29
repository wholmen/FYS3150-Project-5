#include <iostream>
#include <armadillo>
#include <fstream>
#include "functions.h"
#include "pdesolver.h"

using namespace std;
using namespace arma;

void analytical(mat *U);

int main()
{
    int Nx, Nt; double T, dt, h; mat u, u2;

    //---------------------------------------------------------------------------------------
    // Solving 2+1 dimensional diffusion equation with implicit, explicit and analytic method
    //---------------------------------------------------------------------------------------

    initialize(&Nx,&Nt,&T);           // Reading in Nx, Nt and T from user

    h = 1./Nx; dt = T / Nt;           // Spacial steplength h and time step dt
    u = zeros<mat>(Nx,Nx);            // u is a function showing concentration for a given x and y. u(x,y)

    for (int i = 0; i<Nx; i++){       // Setting boundary condition u(0,y) = 1. All other boundarys are 0.
        u(i,0) = 1;
    }
    u2 = u;                           // u2 will be solved by explicit method

    BackwardEuler(&u, Nx, h, dt);     // Solving as a function of h and dt
    ForwardEuler(&u2, Nx, Nt, h, dt); // Solving as a function of h and dt

    ofstream myfile, myfile2;
    myfile.open("Implicit.txt"); myfile2.open("Explicit.txt");
    myfile << u; myfile2 << u2;
    myfile.close(); myfile2.close();

    //---------------------------------------------------------------------------------------------------------
    // Solving 1+1 dimensional diffusion equation using Markov Chain, explicit, implicit and analytical solver.
    //---------------------------------------------------------------------------------------------------------

    int N0; vec X; // N0 is number of particles at x=0. X is the position vector for all particles
    N0 = 1000; X = zeros<vec>(N0); Nt = 1000; dt = 0.01; //All N0 particles start at 0, so X has length N0 with
                                                         //values 0.
    PDESolver D1(Nt,dt);
    D1.RandomWalk(&X);

    myfile.open("RandomWalk_a.txt");
    myfile << X;
    myfile.close();

    return 0;
}

void analytical(mat *U){
    mat u;
    u = *U;


    *U = u;
}

