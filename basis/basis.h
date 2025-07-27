#ifndef basis_H
#define basis_H

#include <iostream>
#include <armadillo>
#include <vector>

class AO{
    arma::vec R;
    arma::vec alpha;
    arma::vec d;
    arma::uvec lmn;
}

class Atom{
    std::vector<AO> AOs;
}

class Molecule{
    std::vector<Atom> Atoms;
    int e;
    int charge;
    int spin_mult;
    Molecule(); //specify parameters
}
