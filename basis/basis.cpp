#include <iostream>
#include <armadillo>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "basis.cpp"

std::map<std::string, int> SYMBOL_MAP {
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

std::map<std::string, std::vector<arma::uvec>> LMN_MAP {
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

Molecule::Molecule(){
    
}



