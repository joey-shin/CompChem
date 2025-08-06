#ifndef input_H
#define input_H

#include <iostream>
#include <armadillo>
#include "../basis/basis.h"

using namespace std;

class Parser{
    private:
        void parse_config(ifstream &filr, map<string, string> &config);
        void parse_molecule(ifstream &file, int &charge, int &spin, vector<string> &atoms, vector<arma::vec> &Rs);
        void parse_Atom(ifstream &file, Atom &atom, arma::vec &R);
        void parse_AO(ifstream &file, const int &K, vector<double> &alphas, vector<double> &ds);
        void parse_SP_AO(ifstream &file, const int &K, vector<double> &alphas, vector<double> &ds_S, vector<double> &ds_P);
    public:
        void parse_config(string &input, map<string, string> &config);
        void parse_molecule(string &input, std::string &basis, Molecule &mol);
        void parse_basis(string &basis_file, vector<string> atoms, vector<arma::vec> Rs, Molecule &mol);
        
};

double FortanDouble(const std::string &str);

#endif
