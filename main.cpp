#include <iostream>
#include <armadillo>

using namespace std;



int main(int argc, char* argv[]){
    
    if (argc < 2){
        printf("Input file not specified!");
        return EXIT_FAILURE;
    }

    string input(argv[1]);

    cout << input;
    return EXIT_SUCCESS;
}
