#include<iostream>
#include<string>
#include"Dispersion.h"
#include<sstream>
#include<fstream>

//Need to initialize the Class Dispersion, first create a constructor

//Constructor
Dispersion::Dispersion() : Natoms_per_monomer(NULL), AtomicSymbols(NULL), 
AtomicNumbers(NULL), spins(NULL), charges(NULL), xyz(NULL), symbols(NULL),
atomic_number(NULL), r(NULL), Symbols(NULL), distance(NULL), other_distance(NULL), xyz_deriv(NULL), xyz_distances(NULL), r2_r4(NULL)
{
}


void Dispersion::Initialize(ifstream &infile) {

    
    CheckInput(infile); 
    NMon = FindNumberOfMonomers(infile);
    CountNumberOfAtoms(infile);
    GetTotalNumberOfAtoms();
    GetCoordinates(infile);
    GetSymbols();
    CalculateAllDistancesAB();
    CalculateDistancesDifferently();
    CoordinateDerivs();
    CoordinateDerivs_2();
    *GetAtomicNumber();
    *GetMultipoleExpectationValue();
   
}

//Destructor
Dispersion::~Dispersion() {

    delete [] spins;
    delete [] charges;
    delete [] Natoms_per_monomer;
    delete [] types;
    delete [] xyz;

    delete [] AtomicSymbols;
    delete [] AtomicNumbers;
    delete [] r2_r4;
    delete [] distance;
    delete [] other_distance;
    delete xyz_deriv;

}    

void Dispersion::CheckInput(ifstream& infile) {
    string line;
    bool start_sxn = false;
    bool end_sxn = false;

    Rewind(infile);
    while (!infile.eof()) { 
        getline(infile, line);
        if (line.substr(0,8) == "molecule") {
            start_sxn = true;
        }
        if (start_sxn == false) {
            cerr << "Calculation failed. The input file is incorrectly formatted." << endl;
            exit(1);
        }
        if (line.substr(0,1) == "}") {
            end_sxn = true;
        }
    }
    infile.clear();
    if (start_sxn=false || end_sxn==false) {
            cerr << "Calculation failed. The input file is incorrectly formatted." << endl;
            exit(1);
    }

}

// This function Counts the number of monomers in the input files
int Dispersion::FindNumberOfMonomers(ifstream& infile) {
    int num_mon;
    string line;
    bool molec_sxn = false;

    // Rewind the file
    Rewind(infile);
    // Count Monomers
    while ( !infile.eof() ) {
        getline(infile,line);

        if (line.substr(0,8) == "molecule") {
            molec_sxn = true;
            num_mon =1;
        }

        if (line.substr(0,2) == "--" && molec_sxn==true) {
            num_mon++;
        }
        
        if (line.substr(0,1) == "}" && molec_sxn==true) {
            break;
        }
    }
    infile.clear();
    return num_mon;
}

//This function counts the number of Atoms in the input file
void Dispersion::CountNumberOfAtoms(ifstream& infile) {
    string line;
    //Rewind the file
    Rewind(infile);
    // Need to create an array that stores the number of atoms per monomer
    Natoms_per_monomer = new int[NMon+1];

    // Create list of total number of monomers and how many atoms per monomer
    int index = 0; // counts the monomers
    int natoms = 0; // counts the atoms within a monomer
    bool molec_sxn = false;
    bool units_declared = false;

    while ( !infile.eof() ) {
        getline(infile,line);
        // Tell the program when it has reached the molecule section
        if (line.substr(0,8)=="molecule") {
            getline(infile,line);
            molec_sxn = true;
        }

        if (line.substr(0,2)=="--" && molec_sxn==true) {
            Natoms_per_monomer[index] = natoms-1;
            natoms = 0;
            getline(infile,line);
            index++;
        }

        if (line.substr(0,14)=="units angstrom" && molec_sxn==true) {
            units_declared = true;
            Natoms_per_monomer[index] = natoms-1;
        }

        if (units_declared == false) {
            Natoms_per_monomer[index] = natoms -1;
        }

        if (line.substr(0,1)=="}" && molec_sxn==true) {
        
            break;
        }
     
        natoms++;
    }
    
}

// Need to know the total number of atoms in the input 
int Dispersion::GetTotalNumberOfAtoms() {
    int Atoms_Total = 0;
    for (int iMon=0;iMon<=NMon-1;iMon++) {
        Atoms_Total += GetNumberOfAtoms(iMon);
    }
    return Atoms_Total;
}

//This is where I get the atom symbols and geometry coordinates
void Dispersion::GetCoordinates(ifstream& infile) {

    string line;
    
    charges = new int[NMon+1];
    spins = new int[NMon+1];
    types = new string[NMon+1];

    int Ntot = GetTotalNumberOfAtoms();
    symbols = new string[Ntot];
    xyz = new double[3*Ntot];

    Rewind(infile);

    int index = 0, atom = 0;
    bool molec_sxn = false;

    while ( !infile.eof() ) {
        getline(infile,line);
        
        if (line.substr(0,8)=="molecule") {
            molec_sxn = true;
            getline(infile,line);

            for (int i=0;i<GetNumberOfAtoms(index);i++) {
            
                infile >> symbols[atom];
                infile >> xyz[3*atom];
                infile >> xyz[3*atom+1];
                infile >> xyz[3*atom+2];
              
                atom++;
            }
            index++;
        }
        if (line.substr(0,2)=="--" && molec_sxn==true) {
            getline(infile,line);
	    
            for (int i=0;i<GetNumberOfAtoms(index);i++) {
                infile >> symbols[atom];
	        infile >> xyz[3*atom];
     		infile >> xyz[3*atom+1];
     		infile >> xyz[3*atom+2];
              
     		atom++;
 	    }
            index++;
        }
        if (line.substr(0,1)=="" && molec_sxn==true) {
        }

        if (line.substr(0,14)=="units angstrom" && molec_sxn==true) {
            break;
        }       
    }
}

double *Dispersion::CalculateAllDistancesAB() {
  
    int Ntot = GetTotalNumberOfAtoms();
    int number_of_distances = 0;
    for (int i=0; i<Ntot; i++) {
        number_of_distances += i;
    }

    double d_x, d_y, d_z;//, r;

    r = new double[number_of_distances];
    int index = -1;

    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++) {
            d_x = xyz[3*i] - xyz[3*j];

            d_y = xyz[3*i+1] - xyz[3*j+1];

            d_z = xyz[3*i+2] - xyz[3*j+2];

            index++;

            r[index] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);

        }
    }
    
    return r;
}

string *Dispersion::GetSymbols() {
    int Ntot = GetTotalNumberOfAtoms();
    Symbols = new string[Ntot];

    for (int i=0; i<Ntot; i++) {
        Symbols[i] = symbols[i];
    }
    return Symbols;
}

double **Dispersion::CalculateDistancesDifferently() {

    int Ntot = GetTotalNumberOfAtoms();
    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;
    distance = new double*[Ntot];
    for(int i=0; i < Ntot; i++) {
        distance[i] = new double[Ntot];
    }

    for (int i=0; i<Ntot;i++ ) {
        for (int j=0; j<Ntot; j++ ) {
            if (i==j) {
                d_x = 0.0;
                d_y = 0.0;
                d_z = 0.0;
            }
            else {
                d_x = xyz[3*i] - xyz[3*j];
                d_y = xyz[3*i+1] - xyz[3*j+1];
                d_z = xyz[3*i+2] - xyz[3*j+2];
            }
            distance[i][j] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
        
            if (i==j) {
                x_deriv = 0.0;
                y_deriv = 0.0;
                z_deriv = 0.0;
            }
            else {
                x_deriv = d_x/distance[i][j];
             
                y_deriv = d_y/distance[i][j];
              
                z_deriv = d_z/distance[i][j];
              
            } 
            
        }
    }
    return distance;
}

double *Dispersion::CoordinateDerivs() {

    int Ntot = GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;
    int index = -1;

    other_distance = new double*[Ntot];
    xyz_deriv = new double[3*combinations];
    for(int i=0; i < Ntot; i++) {
        other_distance[i] = new double[Ntot];
    }

    for (int i=0; i<Ntot;i++ ) {
      
        for (int j=i+1; j<Ntot; j++ ) {
            index+=3;
      
            d_x = xyz[3*i] - xyz[3*j];
          
            d_y = xyz[3*i+1] - xyz[3*j+1];
         
            d_z = xyz[3*i+2] - xyz[3*j+2];
          
            other_distance[i][j] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
         
            xyz_deriv[index-2] = d_x;

            xyz_deriv[index-1] = d_y;

            xyz_deriv[index] = d_z;
            
        }
    }
    for (int i=0; i < 3*combinations; i++) {
   
    }
    return xyz_deriv;
}

double **Dispersion::CoordinateDerivs_2() {

    int Ntot = GetTotalNumberOfAtoms();
    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;

    new_distance = new double*[Ntot];
    for(int i=0; i < Ntot; i++) {
        new_distance[i] = new double[Ntot];
    }
    
    xyz_distances = new double*[Ntot];
    for (int i=0; i < Ntot; i++) {
        xyz_distances[i] = new double[3*Ntot];
    }

    for (int i=0; i<Ntot;i++ ) {
        for (int j=0; j<Ntot; j++ ) {
            if (i==j) {
                d_x = 0.0;
                d_y = 0.0;
                d_z = 0.0;
            }
       
            else {
                d_x = xyz[3*i] - xyz[3*j];
                d_y = xyz[3*i+1] - xyz[3*j+1];
                d_z = xyz[3*i+2] - xyz[3*j+2];
            }
      
            new_distance[i][j] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
            xyz_distances[i][3*j] = d_x;
       
            xyz_distances[i][3*j+1] = d_y;
       
            xyz_distances[i][3*j+2] = d_z;
        
            if (i==j) {
                x_deriv = 0.0;
                y_deriv = 0.0;
                z_deriv = 0.0;
            }
            else {
                x_deriv = d_x/distance[i][j];
             
                y_deriv = d_y/distance[i][j];
             
                z_deriv = d_z/distance[i][j];
             
            } 
            
        }
    }
    return xyz_distances;
}


// C++ does not allow you to return an entire array as an argument to a
// function. To get around this I defined this function as a pointer so
// that the returned argument is a pointer to the array. In this way I can
// access the elements of atomic_number[] from the coord_number class.
int *Dispersion::GetAtomicNumber() {

    int Ntot = GetTotalNumberOfAtoms();
// atomic numbers of elements that we have UCHF coefficients for
// the index of a given element corresponds to its atomic number
    string AtomicSymbols[36] = {
        "XX",
        "H","XX",
        "XX","XX","B","C","N","O","F","Ne", // period 2
        "XX","XX","XX","XX","P","S","Cl","Ar", // period 3
        "XX","XX","XX","XX","XX","XX","XX","XX","XX","XX","XX","XX", // period 4
        "XX","XX","XX","XX","Br"}; // period 4
    
    int i, j, at_num = 0;

    atomic_number = new int[Ntot];
    // The array atomic_number is filled with the atomic numbers of each 
    // atom that appears in the input.
    int counter = 0;
    for (i=0;i<Ntot;i++) {
        for (j=0;j<36;j++) {
            if (symbols[i]==AtomicSymbols[j]) {
                at_num = j;
                atomic_number[i] = at_num;

                counter++;

            }
        }

    }
    if (counter!=Ntot) {
        cerr << "Calculation terminated, the input contains atoms that are not currently supported by MP2D." << endl;
        exit(1);
    }

    // Here the returned argument is a pointer to the array atomic_number[]
    return atomic_number;
}

double *Dispersion::GetMultipoleExpectationValue() {
    int Ntot = GetTotalNumberOfAtoms();

    double MultipoleExpectationValues[36] = { 0.0, 8.0589, 0.0, 0.0, 0.0, 11.8799, 7.8715, 5.5588, 4.7566, 3.8025
,3.1036, 0.0, 0.0, 0.0, 0.0, 9.5361, 8.1652, 6.7463, 5.6004, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.1251};

    int i, j;
    r2_r4 = new double[Ntot]; 

    for (int i=0;i<Ntot;i++) {
    
        r2_r4[i] = MultipoleExpectationValues[atomic_number[i]];
   
    }
    return r2_r4;
} 

