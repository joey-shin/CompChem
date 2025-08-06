#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <armadillo>
#include "basis.h"
#include "../input/input.h"


void AO::printAO() {
    std::cout << "K: " << K << std::endl;
    R.print("R: ");
    alphas.print("alphas: ");
    ds.print("ds: ");
    lmn.print("lmn: ");
    std::cout << "\n";
}


void Atom::addAO(AO &ao) {
    AOs.push_back(ao);
}


void Atom::printAtom() {
    std::cout << "Symbol: " << symbol << "\nAtomic Number: " << z << "\n" << std::endl;
    for (auto ao : AOs) {
        ao.printAO();
    }
}


Molecule::Molecule(std::string &input, std::string &basis, Parser &parser){
    std::vector<std::string> atoms;
    std::vector<arma::vec> Rs;

    parser.parse_molecule(input, charge, spin, atoms, Rs); // parse %molecule input
    std::string basis_file = "basis/basis/" + basis + ".bas";
    parser.parse_basis(basis_file, atoms, Rs, *this);

}


void Molecule::printMolecule() {
    std::cout << "Charge: " << charge << "\nSpin: " << spin << "\nNumber of Atoms: " << Atoms.size() << "\n" <<std::endl;;
    for (auto atom : Atoms) {
        atom.printAtom();
    }
}


void Molecule::addAtom(Atom &atom) {
    Atoms.push_back(atom);
}


