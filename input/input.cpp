#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <armadillo>
#include "input.h"

using namespace std;

void Parser::parse_config(ifstream &file, map<string, string> &config) { // parse %config settings
    string line;
    while (getline(file, line)) {
        if (line.find("%end") != string::npos) {
            break;
        }

        string element, option;
        istringstream iss(line);
        if (!(iss >> element >> option)){
            throw invalid_argument("invalid %config format");
        }

        config.insert({element, option});
    }
}


void Parser::parse_molecule(ifstream &file, int &charge, int &spin, vector<string> &atoms, vector<arma::vec> &Rs) { // parse %molecule input
    string line;
    getline(file, line);
    istringstream iss(line);

    if (!(iss >> charge >> spin)){
        throw invalid_argument("Invalid charge/spin format in %Molecule");
    }

    while (getline(file, line)) {
        istringstream iss(line);
        if (line.find("%end") != string::npos) {
            break;
        }

        string atom;
        arma::vec R(3);
        if (!(iss >> atom >> R[0] >> R[1] >> R[2])){
            throw invalid_argument("Invalid Atom/Coord format in %molecule");
        }
        atoms.push_back(atom);
        Rs.push_back(R);
    }
}



void Parser::parse_config(string &input, map<string, string> &config) { // public invocation of parse_config
    ifstream file(input);
    if (!file.is_open()) {
        throw invalid_argument("Could not open input file: " + input);
    }

    string line;
    while (getline(file, line)) {
        if (line.find("%config") != string::npos) {
            parse_config(file, config);
        }
    }
}


void Parser::parse_molecule(string &input, std::string &basis, Molecule &mol) { // public invocation of parse_molecule
    ifstream file(input);
    if (!file.is_open()) {
        throw invalid_argument("Could not open input file: " + input);
    }

    string line;
    vector<string> atoms;
    vector<arma::vec> Rs;
    while (getline(file, line)) {
        if (line.find("%molecule") != string::npos) {
            parse_molecule(file, mol.charge, mol.spin, atoms, Rs);
        }
    }
    string basis_file = "basis/basis/" + basis + ".bas";
    parse_basis(basis_file, atoms, Rs, mol);
}


const std::map<std::string, int> SYMBOL_MAP {
    {"X", 0}, {"H", 1}, {"He", 2}, // "X": ghost atom
    {"Li", 3}, {"Be", 4}, 
    {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, 
    {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
    {"K", 19}, {"Ca", 20}, 
    {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
    {"Rb", 37}, {"Sr", 38}, 
    {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
    {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54},
    {"Cs", 55}, {"Ba", 56}, 
    {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71},
    {"Hf", 72}, {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}
};


const std::map<std::string, std::vector<arma::uvec>> LMN_MAP {
    {"S", {arma::uvec({0, 0, 0})}},
    {"P", {arma::uvec({1, 0, 0}), 
           arma::uvec({0, 1, 0}), 
           arma::uvec({0, 0, 1})}},
    {"D", {arma::uvec({2, 0, 0}), 
           arma::uvec({1, 1, 0}), 
           arma::uvec({1, 0, 1}),
           arma::uvec({0, 2, 0}),
           arma::uvec({0, 1, 1}),
           arma::uvec({0, 0, 2})}},
    
    {"F", {arma::uvec({3, 0, 0}),
           arma::uvec({2, 1, 0}),
           arma::uvec({2, 0, 1}),
           arma::uvec({1, 2, 0}),
           arma::uvec({1, 1, 1}),
           arma::uvec({1, 0, 2}),
           arma::uvec({0, 3, 0}),
           arma::uvec({0, 2, 1}),
           arma::uvec({0, 1, 2}),
           arma::uvec({0, 0, 3})}}
};


double FortanDouble(const string &str){
    std::string fix = str;
    std::replace(fix.begin(), fix.end(), 'D', 'E');
    return std::stod(fix);
}


bool isint(const string &s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}


bool isdouble(const string &s) {
    std::istringstream iss(s);
    double val;
    char nonval;
    return (iss >> val) && !(iss >> nonval);
}


bool isShell(const string &s) {
    const std::vector<std::string> shellLabels = {"S", "P", "D", "F", "G", "SP"};
    return std::find(shellLabels.begin(), shellLabels.end(), s) != shellLabels.end();
}



void Parser::parse_SP_AO(ifstream &file, const int &K, vector<double> &alphas, vector<double> &ds_S, vector<double> &ds_P) {
    string line;
    string first, second, third;

    for (int i = 0; i < K; i++) {
        getline(file, line);
        istringstream iss(line);
        if (!(iss >> first >> second >> third)) {
            throw invalid_argument("Invalid AO format in basis file");
        }
        alphas.push_back(FortanDouble(first));
        ds_S.push_back(FortanDouble(second));
        ds_P.push_back(FortanDouble(third));
    }
 
}


void Parser::parse_AO(ifstream &file, const int &K, vector<double> &alphas, vector<double> &ds) { // parse exponents and contraction coefficients
    string line;
    string first, second;
    
    for (int i = 0; i < K; i++){
        getline(file, line);
        istringstream iss(line);
        if (!(iss >> first >> second)) {
            throw invalid_argument("Invalid AO format in basis file");
        }
        alphas.push_back(FortanDouble(first));
        ds.push_back(FortanDouble(second));
    }
}



void Parser::parse_Atom(ifstream &file, Atom &atom, arma::vec &R) {
    string line;
    string first, second, third;
    string shell;
    int K; //number of pGTOs
    // double scaling; // scaling factor

    while (getline(file, line)) {
        istringstream iss(line);
        iss >> first >> second >> third;
        
        if (first == "****") {
            break;
        } else if (isShell(first) && isint(second) && isdouble(third)) {
            K = stoi(second);
            // scaling = FortanDouble(third);
            if (first == "SP") { // create AO for S and P seperately for SP shell
                vector<double> alphas, ds_S, ds_P;
                parse_SP_AO(file, K, alphas, ds_S, ds_P);

                AO ao_S(K, R, arma::vec(alphas), arma::vec(ds_S), arma::uvec({0, 0, 0}));
                atom.addAO(ao_S);

                const auto& lmns = LMN_MAP.at("P");
                for (const arma::uvec &lmn : lmns) {
                    AO ao_P(K, R, arma::vec(alphas), arma::vec(ds_P), lmn);
                    atom.addAO(ao_P);
                }

            } else if (LMN_MAP.find(first) != LMN_MAP.end()) { // rest of the shell types
                const auto& lmns = LMN_MAP.at(first);
                vector<double> alphas, ds;
                parse_AO(file, K, alphas, ds);

                for (const arma::uvec &lmn : lmns) {
                    AO ao(K, R, arma::vec(alphas), arma::vec(ds), lmn);
                    atom.addAO(ao);
                }
            } else {
                throw invalid_argument("Invalid AO format in basis file");
            }
        } 
    }
}


void Parser::parse_basis(string &basis_file, vector<string> atoms, vector<arma::vec> Rs, Molecule &mol) {
    std::vector<int> zs;
    for (const std::string &symbol : atoms) { // construct vector of atomic numbers
        auto z = SYMBOL_MAP.find(symbol);
        if (z == SYMBOL_MAP.end()) {
            throw std::invalid_argument("Unknown atom symbol in %molecule: " + symbol);
        }
        zs.push_back(z->second);
    }

    //sum up to calculate number of electrons

    string line;
    for (int i = 0; i < atoms.size(); i++) { // loop through each atom
        ifstream file(basis_file);
        if (!file.is_open()){
            throw invalid_argument("Could not open basis file: " + basis_file);
        }
        while (getline(file, line)){
            istringstream iss(line);
            string first, second, third;
            iss >> first >> second >> third;

            if (first == atoms[i] && isint(second) && third == "") { // search if atom index i from atoms matches symbol in basis
                // int ECP = stoi(second); // store effective core potential
                Atom atom;
                atom.symbol = first;
                atom.z = SYMBOL_MAP.at(first);
                parse_Atom(file, atom, Rs[i]);
                mol.addAtom(atom);
                break;
            }
        }
    }
}


