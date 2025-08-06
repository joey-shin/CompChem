#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <armadillo>
#include <vector>


struct AO { // Each AO (cGTO) contains Contraction Coeff (d) and Gauss Exp (alpha) of pGTOs
    AO(const int &K, const arma::vec &R, const arma::vec &alphas, const arma::vec &ds, const arma::uvec &lmn)
        : K(K), R(R), alphas(alphas), ds(ds), lmn(lmn) {}
    int K; // number og pGTOs
    arma::vec R;
    arma::vec alphas; // pGTO exponent
    arma::vec ds; // contraction coeff
    arma::uvec lmn;
    void printAO();
};

struct Atom {
    std::vector<AO> AOs;
    std::string symbol; // atomic symbol
    int z; // atomic number
    void addAO(AO &ao);
    void printAtom();
};

struct Molecule {
    std::vector<Atom> Atoms;
    int e;
    int charge;
    int spin;
    void addAtom(Atom &atom);
    void printMolecule();
};

#endif
