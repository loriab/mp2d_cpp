#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>

//using std::ifstream;
//using std::string;
//#include<Coordination_Number>

using namespace std;

// The dispersion class will be the main body of code. This section will parse
// the input to obtain the necessary parameters from the geometry.
// The parameters will be passed to the other modules to calculate CN, Damping
// function, etc. Then the results will be passed back here to calculate E6 and E8

class Dispersion {
 private:
    Dispersion();
    ~Dispersion();

    int NMon; // number of monomers
    // Pointers for arrays to store info about each monomer. 
    int *Natoms_per_monomer, *spins, *charges;
    string *types;
    string *AtomicSymbols;
    string *symbols;
    int *AtomicNumbers;
    int *atomic_number;
    double *r2_r4;
    // I'm making the array for coordinates global, but it is local in HMBI
    double *xyz;
//    double *r;
   
    //define some variables that will be used to calculate the dispersion energy
    double C6;
    double C8;
    double E6;
    double E8;
    double S6;
    double S8;

    // I also need to define the function that will be used to calculate the dispersion energy
 public:
    // I need to instantiate the single instance of this class
    static Dispersion& dispersion() {
        static Dispersion thedispersion;
        return thedispersion;
    //    return 0;
    }
    void Initialize(ifstream& infile);
    // I need to know how many molecules are in the input file
    int FindNumberOfMonomers(ifstream& infile);
    // I need to know how many atoms are in the input file
    void CountNumberOfAtoms(ifstream& infile);
    // Need to be able to get the number of atoms in each monomer for various other functions
    int GetNumberOfAtoms(int iMon) {return Natoms_per_monomer[iMon];};
    // Also need to get the total number of atoms in the whole input.
    int GetTotalNumberOfAtoms(); // {return Natoms_per_monomer[0];};

    //Routines for getting coordinates, atom symbols and other stuff from the input file
    void SetAtomicSymbols(); // I don't think I need this one
    void GetCoordinates(ifstream& infile);
    string *GetSymbols();
//    void CalculateAllDistancesAB();
    double * CalculateAllDistancesAB();
    double **CalculateDistancesDifferently();
// Original CoordinateDerivs Function before 2-D distances array
    double *CoordinateDerivs();
    double **CoordinateDerivs_2();
// New CoordinateDerivs function to work with Under Development 2
//    double **CoordinateDerivs();
    int *GetAtomicNumber();
    double *GetMultipoleExpectationValue();
    int GetCovalentRadii();
    double GetC6();

    double *r;
    string *Symbols;
    double **distance;
    double **other_distance;
    double **new_distance;
    // Original xyz_deriv before 2-D distances array
    double *xyz_deriv;
    double **xyz_distances;
//New xyz_deriv to work with Under Development 2
  //  **xyz_deriv;

    double E6_Dispersion(); 
    // This function will calculate the R^6 dependent Dispersion Energy Contribution.

    double E8_Dispersion(); 
    // This function will calculate the R^8 dependent Dispersion Energy Contribution

    //int GetNumberOfMonomers() {return NMon;};
    // Not sure if all of this needs to be here in rewind
    void Rewind(ifstream& infile) {
        infile.clear();
        infile.seekg(0,ios::beg);
    }
};
