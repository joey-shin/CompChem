#include <iostream>
#include <armadillo>
#include <math.h>

using namespace std;

double E(int i, int j, int t, double Qx, double a, double b){
    // a and b: guassian exponents
    // i and j: angular momentum number for gaussians a and b
    // Qx: distance between gaussians a and b
    
    double p = a + b; 
    double q = a * b / (a + b); 
    if (t < 0 || t > i + j){
        return 0.0;
    }
    if (i == 0 && j == 0 && t ==0){
        return exp(-q * Qx * Qx);
    }
    
    if (j == 0){ // terminate j, increment i
        return (1 / (2 * p)) * E(i-1, j, t-1, Qx, a, b) -
               (q * Qx / a) * E(i-1, j, t, Qx, a, b) +
               (t + 1) * E(i-1, j, t+1, Qx, a, b);
    }
    else{ // increment j if j isn't terminus
        return (1 / (2 * p)) * E(i, j-1, t-1, Qx, a, b) + 
               (q * Qx / a) * E(i, j-1, t, Qx, a, b) +
               (t + 1) * E(i, j-1, t+1, Qx, a, b);
    }
}

double S_uv(arma::uvec lmn_a, arma::uvec lmn_b, double a, double b, arma::vec A, arma::vec B){ //pGTO overlap (no multipole moments)
    // lmn_a and lmn_b: vector of angular momentum numbers
    // a and b: gaussian exponents
    // i and j: angular moment
    // A and B: 3d coordinates for gaussian centers

    double Sx = E(lmn_a(0), lmn_b(0), 0, A(0) - B(0), a, b);
    double Sy = E(lmn_a(1), lmn_b(1), 0, A(1) - B(1), a, b);
    double Sz = E(lmn_a(2), lmn_b(2), 0, A(2) - B(2), a, b);
    
    return Sx * Sy * Sz * pow((M_PI / (a + b)), 1.5); 
}


