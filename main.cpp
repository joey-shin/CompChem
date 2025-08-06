#include <iostream>
#include <armadillo>
#include <stdexcept>

#include "input/input.h"
#include "basis/basis.h"


int main(int argc, char* argv[]){
    
    printf("\n\nv0.2\n");
    printf("A Confusing Computational Chemistry Program made by a Confused Undergrad\n");
    printf("\"i'm confused send help\"\n\n");

    if (argc < 2){
        std::printf("Input file not specified!");
        return EXIT_FAILURE;
    }

    std::string input(argv[1]); // get input file 
    std::cout << "Input file accepted: " << input << std::endl;

    try {
        std::printf("\n***\nLoading input module...\n\n");
        Parser parser; 
        std::map<std::string, std::string> config;
        parser.parse_config(input, config); // parse configuration settings

        std::printf("\n***\nLoading basis module...\n\n");
        
        Molecule mol;
        parser.parse_molecule(input, config["basis:"], mol); // take in molecule object and directly construct Molecule
        // mol.printMolecule();

    } catch (std::invalid_argument &e) {
        std::cerr << "\n\nCaught Exception std::invalid_argument: " << e.what() << endl;
        return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;
}
