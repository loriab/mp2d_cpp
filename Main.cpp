#include<iostream>
#include<string>
#include"Dispersion.h"
#include"Coord_Num.h"
#include<sstream>
#include<fstream>
#include <cstring>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

int main(int argc, char ** argv) {
    // This opens the input file
    const int maxlength = 50;
    char filename[maxlength];
    if (!argv[1]) {
        cerr << "No input file provided" << endl;
        return 1;
    }
    strcpy(filename, argv[1]);
    
    ifstream infile;
    infile.open(filename);
//HMBI uses this line but originally I did not so comment assert() out if you wish to return to the true original
    assert(infile.is_open());
//    infile.close();

    // I have to initialize the Dispersion object
    //Dispersion::Dispersion().Initialize(infile);
    Dispersion::dispersion().Initialize(infile);
    //int NMon = Dispersion::dispersion().GetNumberOfMonomers();
    // I need to initialize the Coord_Num object
    Coord_Num::coord_num().Initialize(infile);

}

