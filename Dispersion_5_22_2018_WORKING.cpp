#include<iostream>
#include<string>
#include"Dispersion.h"
#include<sstream>
#include<fstream>

//using namespace std;
//This will be the meat of the MP2D Calculation
//This is where I will parse the input to grab the number of atoms, atom type,
//number of monomers, and coordinates.

//Need to initialize the Class Dispersion, first create a constructor

//Constructor
Dispersion::Dispersion() : Natoms_per_monomer(NULL), AtomicSymbols(NULL), 
AtomicNumbers(NULL), spins(NULL), charges(NULL), xyz(NULL), symbols(NULL),
atomic_number(NULL), r(NULL), Symbols(NULL), distance(NULL), other_distance(NULL), xyz_deriv(NULL), xyz_distances(NULL), r2_r4(NULL)
{
}


void Dispersion::Initialize(ifstream &infile) {

    
    NMon = FindNumberOfMonomers(infile);
    CountNumberOfAtoms(infile);
    GetTotalNumberOfAtoms();
    GetCoordinates(infile);
    GetSymbols();
    CalculateAllDistancesAB();
    CalculateDistancesDifferently();
    CoordinateDerivs();
    CoordinateDerivs_2();
//    SetAtomicSymbols(); // I don't think I'm going to need this one
    *GetAtomicNumber();
    *GetMultipoleExpectationValue();
//    GetCovalentRadii();
   
}

//Destructor
Dispersion::~Dispersion() {

    delete [] spins;
    delete [] charges;
    delete [] Natoms_per_monomer;
    delete [] types;
    delete [] xyz;
//    delete [] symbols;

    delete [] AtomicSymbols;
    delete [] AtomicNumbers;
    delete [] r2_r4;
    delete [] distance;
    delete [] other_distance;
    delete xyz_deriv;

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
 //   cout << "FindNumberOfMonomers_Output is: " <<  num_mon << endl;
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
      //      cout << line << endl;
            index++;
     //       cout << "the index is: " << index << endl;
        }

        if (line.substr(0,14)=="units angstrom" && molec_sxn==true) {
            units_declared = true;
            Natoms_per_monomer[index] = natoms-1;
        }

        if (units_declared == false) {
            Natoms_per_monomer[index] = natoms -1;
        }

        if (line.substr(0,1)=="}" && molec_sxn==true) {
         //   Natoms_per_monomer[index] = natoms-1;
            break;
        }
     //   else { getline(infile, line);}
        natoms++;
    }
    for (int i=0;i<=NMon-1;i++) {
       // cout << "Now NMon is: " << NMon << endl;
       // cout << Natoms_per_monomer[i] << endl;
  //      cout << "Monomer " << i << " has " << Natoms_per_monomer[i] << " atoms " << endl; 
    }
    // This is just a test to print out the Number of Atoms in each Monomer
    for (int iMon=0;iMon<=NMon-1;iMon++) {
        int N = GetNumberOfAtoms(iMon); 
    //    cout << N << endl;
    }
//    int G = GetTotalNumberOfAtoms();
//    cout << G << endl;
}
// Need to know the total number of atoms in the input 
int Dispersion::GetTotalNumberOfAtoms() {
    int Atoms_Total = 0;
    for (int iMon=0;iMon<=NMon-1;iMon++) {
        Atoms_Total += GetNumberOfAtoms(iMon);
    }
  //  cout << Atoms_Total << endl;
    return Atoms_Total;
}


// Here is where I'm going to try and get the atomic coordinates
// when you left to go judge posters this function appears nowhere else but here and
// it is completely experimental. There is no reason to think any of this will work.
/*void Dispersion::GetCoordinates(ifstream& infile) {
    string line;

    Rewind(infile);
    cout << "is there anybody out there" << endl;
    string a;
    double b;
    double c;
    double d;

    string element;
    double x_coordinate;
    double y_coordinate;
    double z_coordinate;

    while (!infile.eof() ) {
        getline(infile,line);
        if (line.substr(0,8)=="molecule") {
            getline(infile,line);
            cout << "hey hey" << endl;
        }
        if (infile>>a>>b>>c>>d)  {
            element = a;
//            cout << element << endl;
            x_coordinate = b;
//            cout << x_coordinate << endl;
            y_coordinate = c;
//            cout << y_coordinate << endl;
            z_coordinate = d;
//            cout << z_coordinate << endl;
//            printf ("%s%5f%5f%5f", element.cstr(), x_coordinate, y_coordinate, z_coordinate);
            cout << element << "     " << x_coordinate << "     " << y_coordinate << "     " << z_coordinate << "     " << endl;
        }
        
        else {
            getline(infile,line);
        }
       /* if (line.substr(0,2)=="--") {
            getline(infile,line);
        }
        if (line.substr(0,1)=="}") {
            break;
        }
    }
}*/

//This is my second attempt to get the atom symbols and geometry coordinates
// This Works!
void Dispersion::GetCoordinates(ifstream& infile) {

    string line;
    
    charges = new int[NMon+1];
    spins = new int[NMon+1];
    types = new string[NMon+1];

    int Ntot = GetTotalNumberOfAtoms();
//    string symbols[Ntot];
    symbols = new string[Ntot];
//    double *xyz;
    xyz = new double[3*Ntot];

    Rewind(infile);

    int index = 0, atom = 0;
    bool molec_sxn = false;

    while ( !infile.eof() ) {
        getline(infile,line);
        
        if (line.substr(0,8)=="molecule") {
            molec_sxn = true;
            getline(infile,line);
            //index++;
//            cout << "index:" << index << endl;
//            cout << line << endl;
//            cout << GetNumberOfAtoms(0) << endl;
            for (int i=0;i<GetNumberOfAtoms(index);i++) {
            
                infile >> symbols[atom];
		//cout << symbols[atom] << endl;
                infile >> xyz[3*atom];
                infile >> xyz[3*atom+1];
                infile >> xyz[3*atom+2];
              //  cout << symbols[atom] << "  " << xyz[3*atom] << "  " << xyz[3*atom+1] << "  " << xyz[3*atom+2] << endl;
                atom++;
            }
            index++;
        }
        if (line.substr(0,2)=="--" && molec_sxn==true) {
            getline(infile,line);
	    //index++;
            for (int i=0;i<GetNumberOfAtoms(index);i++) {
                infile >> symbols[atom];
//     		cout << symbols[atom] << endl;
	        infile >> xyz[3*atom];
     		infile >> xyz[3*atom+1];
     		infile >> xyz[3*atom+2];
              //  cout << symbols[atom] << "  " << xyz[3*atom] << "  " << xyz[3*atom+1] << "  " << xyz[3*atom+2] << endl;
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
    // Just a test print to make sure I can access these arrays
   /* for (int i=0;i<Ntot;i++) {
        cout << symbols[i] << endl;
        cout << xyz[i] << endl; 
    } */
}

//void Dispersion::CalculateAllDistancesAB() {
double *Dispersion::CalculateAllDistancesAB() {
    // Can I calculate the distances between atom pairs AB?
// Yes I Can!
/*    double x_diff;
    x_diff = xyz[1] - xyz[4];
    cout << x_diff << endl;*/

    int Ntot = GetTotalNumberOfAtoms();
    int number_of_distances = 0;
    for (int i=0; i<Ntot; i++) {
        number_of_distances += i;
    }
//    cout << "the number of distances is" << number_of_distances << endl;
    double d_x, d_y, d_z;//, r;
//    double r[number_of_distances];
    r = new double[number_of_distances];
    int index = -1;
//    cout << Ntot << endl;
    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++) {
            d_x = xyz[3*i] - xyz[3*j];
//	    cout << d_x << endl;
            d_y = xyz[3*i+1] - xyz[3*j+1];
//            cout << d_y << endl;
            d_z = xyz[3*i+2] - xyz[3*j+2];
//            cout << d_x << "  " << d_y << "  " << d_z << endl;
            index++;

            r[index] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
//            cout << r[index] << endl;
        }
    }
    for (int i=0;i<number_of_distances;i++) {
  //      cout << r[i] << endl;
    }
    return r;
}

string *Dispersion::GetSymbols() {
    int Ntot = GetTotalNumberOfAtoms();
    Symbols = new string[Ntot];

    for (int i=0; i<Ntot; i++) {
        Symbols[i] = symbols[i];
//        cout << Symbols[i] << endl;
    }
    return Symbols;
}

double **Dispersion::CalculateDistancesDifferently() {
//    cout << "is there anybody out there" << endl;
    int Ntot = GetTotalNumberOfAtoms();
    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;
//    double distance[Ntot][Ntot];
//    distance = new double[Ntot][Ntot];
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
        //    cout << "the new distance is " << distance[i][j] << endl;
            if (i==j) {
                x_deriv = 0.0;
                y_deriv = 0.0;
                z_deriv = 0.0;
            }
            else {
                x_deriv = d_x/distance[i][j];
              //  cout << "x_deriv " << x_deriv << endl;
                y_deriv = d_y/distance[i][j];
              //  cout << "y_deriv " << y_deriv << endl;
                z_deriv = d_z/distance[i][j];
              //  cout << "z_deriv " << z_deriv << endl;
            } 
            
        }
    }
    return distance;
}

double *Dispersion::CoordinateDerivs() {
//    cout << "is there anybody out there" << endl;

    int Ntot = GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;
    int index = -1;
//    double distance[Ntot][Ntot];
//    distance = new double[Ntot][Ntot];
    other_distance = new double*[Ntot];
    xyz_deriv = new double[3*combinations];
    for(int i=0; i < Ntot; i++) {
        other_distance[i] = new double[Ntot];
    }

    for (int i=0; i<Ntot;i++ ) {
       // for (int j=0; j<Ntot; j++ ) {
        for (int j=i+1; j<Ntot; j++ ) {
            index+=3;
       //     cout << index << endl;
          //  if (i==j) {
          //      d_x = 0.0;
          //      d_y = 0.0;
          //      d_z = 0.0;
          //  }
          //  else {
            d_x = xyz[3*i] - xyz[3*j];
          //  cout << "d_x " << d_x << endl;
            d_y = xyz[3*i+1] - xyz[3*j+1];
          //  cout << "d_y " << d_y << endl;
            d_z = xyz[3*i+2] - xyz[3*j+2];
           // cout << "d_z " << d_z << endl;
         //   }
            other_distance[i][j] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
         //   cout << "the new distance is " << other_distance[i][j] << endl;
            //if (i==j) {
            //    x_deriv = 0.0;
             //   cout << "x_deriv " << x_deriv << endl;
            //    y_deriv = 0.0;
            //    cout << "y_deriv " << y_deriv << endl;
            //    z_deriv = 0.0;
            //    cout << "z_deriv " << z_deriv << endl;
           // }
          //  else {
//            x_deriv = d_x/other_distance[i][j];
       //     x_deriv = xyz[3*i]/other_distance[i][j];
       //     cout << "x_deriv " << x_deriv << endl;
       //     cout << other_distance[i][j] << endl;
//            y_deriv = d_y/other_distance[i][j];
       //     y_deriv = xyz[3*i+1]/other_distance[i][j];
       //     cout << "y_deriv " << y_deriv << endl;
//            z_deriv = d_z/other_distance[i][j];
       //     z_deriv = xyz[3*i+2]/other_distance[i][j];
       //     cout << "z_deriv " << z_deriv << endl;
//            xyz_deriv[index-2] = x_deriv;
            xyz_deriv[index-2] = d_x;
//            xyz_deriv[index-1] = y_deriv;
            xyz_deriv[index-1] = d_y;
//            xyz_deriv[index] = z_deriv;
            xyz_deriv[index] = d_z;
          //  } 
            
        }
    }
    for (int i=0; i < 3*combinations; i++) {
    //    cout << xyz_deriv[i] << endl;
    }
    return xyz_deriv;
}

double **Dispersion::CoordinateDerivs_2() {
//    cout << "is there anybody out there" << endl;
    int Ntot = GetTotalNumberOfAtoms();
    double d_x, d_y, d_z;
    double x_deriv, y_deriv, z_deriv;
//    double distance[Ntot][Ntot];
//    distance = new double[Ntot][Ntot];
    new_distance = new double*[Ntot];
    for(int i=0; i < Ntot; i++) {
        new_distance[i] = new double[Ntot];
    }
    //   xyz_distances = new double[3*combinations];
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
         //   else if (i<j) {
            else {
                d_x = xyz[3*i] - xyz[3*j];
                d_y = xyz[3*i+1] - xyz[3*j+1];
                d_z = xyz[3*i+2] - xyz[3*j+2];
            }
       //     else if (i>j) {
         //       d_x = xyz[3*j] - xyz[3*i];
           //     d_y = xyz[3*j+1] - xyz[3*i+1];
         //       d_z = xyz[3*j+2] - xyz[3*i+2];
      //      }
      //      else {}
            new_distance[i][j] = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
            xyz_distances[i][3*j] = d_x;
        //    cout << xyz_distances[i][3*j] << endl;
            xyz_distances[i][3*j+1] = d_y;
        //    cout << xyz_distances[i][3*j+1] << endl;
            xyz_distances[i][3*j+2] = d_z;
         //   cout << xyz_distances[i][3*j+2] << endl;
         //   cout << "the new distance is " << new_distance[i][j] << endl;
       //     cout << d_x << endl;
       //     cout << d_y << endl;
       //     cout << d_z << endl;
            if (i==j) {
                x_deriv = 0.0;
                y_deriv = 0.0;
                z_deriv = 0.0;
            }
            else {
                x_deriv = d_x/distance[i][j];
              //  cout << "x_deriv " << x_deriv << endl;
                y_deriv = d_y/distance[i][j];
              //  cout << "y_deriv " << y_deriv << endl;
                z_deriv = d_z/distance[i][j];
              //  cout << "z_deriv " << z_deriv << endl;
            } 
            
        }
    }
    return xyz_distances;
}




// Need to create arrays that store atomic numbers and atomic symbols
//Need to implement an error checking mechanism for input files that contain atoms for which
//we do not have UCHF Coefficients
/*void Dispersion::SetAtomicSymbols() {
//    AtomicSymbols = new string[GetTotalNumberOfAtoms()];
//    AtomicNumbers = new int[GetTotalNumberOfAtoms()];
    cout << "SetAtomicSymbols" << endl;

}*/

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
//    int atomic_number[Ntot];
    atomic_number = new int[Ntot];
    // The array atomic_number is filled with the atomic numbers of each 
    // atom that appears in the input.
    for (i=0;i<Ntot;i++) {
        for (j=0;j<36;j++) {
            if (symbols[i]==AtomicSymbols[j]) {
                at_num = j;
                atomic_number[i] = at_num;
//		cout << atomic_number[i] << endl;
//		cout << "atom  " << symbols[i] << " has atomic number: " << at_num << endl;
            }
        }
//    cout << atomic_number[i] << endl;
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


 /*   for (int i=0;i<Ntot;i++) {
        for (int j=0;j<Ntot;j++) {
            cout << MultipoleExpectationValues[atomic_number[j]] << endl;
        }
    } */
    for (int i=0;i<Ntot;i++) {
    //    cout << MultipoleExpectationValues[atomic_number[i]] << endl;
        r2_r4[i] = MultipoleExpectationValues[atomic_number[i]];
    //    cout << r2_r4[i] << endl;
    }
    return r2_r4;
} 

/*int Dispersion::GetCovalentRadii() {
    // Covalent radii for H,B,C,N,O,F,Ne,P,S,Cl,Ar, and Br
    // Each covalent radii is placed at the index corresponding to the element's
    // atomic number
    double CovalentRadii[36] =  {
        0.0,
        0.32,0.0,
        0.0,0.0,0.77,0.75,0.71,0.63,0.64,0.67, // period 2
        0.0,0.0,0.0,0.0,1.1,1.02,0.99,0.96, // period 3
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // period 4
        0.0,0.0,0.0,0.0,1.14}; // period 4

    int Ntot = GetTotalNumberOfAtoms();

    for (int i=0; i<Ntot; i++ ) {
//        cout << atomic_number[i] << endl;
        cout << CovalentRadii[atomic_number[i]] << endl;
    }
} */

/*int Dispersion::GetC6() {
    
    // This is where I get the precomputed Dispersin Coefficients.
    double C6_Matrix[2975] = {
1,	1,	1,	1,	2.402,
1,	1,	0,	1,	2.9731,
1,	1,	0,	0,	3.6967,
1,	5,	1,	4,	8.8628,
1,	5,	1,	3,	8.4651,
1,	5,	1,	2,	9.8734,
1,	5,	1,	1,	10.9649,
1,	5,	1,	0,	14.0332,
1,	5,	0,	4,	11.0123,
1,	5,	0,	3,	10.5109,
1,	5,	0,	2,	12.3093,
1,	5,	0,	1,	13.7008,
1,	5,	0,	0,	17.6137,
1,	6,	1,	4,	7.5675,
1,	6,	1,	3,	9.1573,
1,	6,	1,	2,	9.9431,
1,	6,	1,	1,	9.4839,
1,	6,	1,	0,	9.6245,
1,	6,	0,	4,	9.3414,
1,	6,	0,	3,	11.3732,
1,	6,	0,	2,	12.3661,
1,	6,	0,	1,	11.7733,
1,	6,	0,	0,	11.9483,
1,	7,	1,	3,	6.724,
1,	7,	1,	2,	7.7094,
1,	7,	1,	1,	6.7637,
1,	7,	1,	0,	6.9518,
1,	7,	0,	3,	8.2699,
1,	7,	0,	2,	9.5104,
1,	7,	0,	1,	8.3178,
1,	7,	0,	0,	8.5545,
1,	8,	1,	2,	5.0379,
1,	8,	1,	1,	5.164,
1,	8,	1,	0,	5.3584,
1,	8,	0,	2,	6.1416,
1,	8,	0,	1,	6.3024,
1,	8,	0,	0,	6.5496,
1,	9,	1,	1,	3.7454,
1,	9,	1,	0,	4.1786,
1,	9,	0,	1,	4.5279,
1,	9,	0,	0,	5.0775,
1,	10,	1,	0,	3.3122,
1,	10,	0,	0,	4.0045,
1,	15,	1,	3,	21.3001,
1,	15,	1,	2,	25.167,
1,	15,	1,	1,	21.0517,
1,	15,	1,	0,	21.5232,
1,	15,	0,	3,	26.6507,
1,	15,	0,	2,	31.6241,
1,	15,	0,	1,	26.3372,
1,	15,	0,	0,	26.9383,
1,	16,	1,	2,	20.0701,
1,	16,	1,	1,	19.1618,
1,	16,	1,	0,	18.3441,
1,	16,	0,	2,	25.024,
1,	16,	0,	1,	23.871,
1,	16,	0,	0,	22.831,
1,	17,	1,	1,	16.7804,
1,	17,	1,	0,	15.6893,
1,	17,	0,	1,	20.8029,
1,	17,	0,	0,	19.4236,
1,	18,	1,	0,	13.4894,
1,	18,	0,	0,	16.6191,
1,	35,	1,	1,	23.2356,
1,	35,	1,	0,	21.4927,
1,	35,	0,	1,	28.9071,
1,	35,	0,	0,	26.6989,
5,	5,	4,	4,	32.81,
5,	5,	4,	3,	31.319,
5,	5,	4,	2,	36.6541,
5,	5,	4,	1,	40.7833,
5,	5,	4,	0,	52.394,
5,	5,	3,	3,	29.8992,
5,	5,	3,	2,	34.9715,
5,	5,	3,	1,	38.8982,
5,	5,	3,	0,	49.9389,
5,	5,	2,	2,	41.0595,
5,	5,	2,	1,	45.7643,
5,	5,	2,	0,	58.9957,
5,	5,	1,	1,	51.0657,
5,	5,	1,	0,	65.977,
5,	5,	0,	0,	85.6176,
5,	6,	4,	4,	27.8585,
5,	6,	4,	3,	33.8852,
5,	6,	4,	2,	36.8364,
5,	6,	4,	1,	35.0808,
5,	6,	4,	0,	35.6024,
5,	6,	3,	4,	26.6205,
5,	6,	3,	3,	32.3498,
5,	6,	3,	2,	35.1593,
5,	6,	3,	1,	33.4927,
5,	6,	3,	0,	33.9901,
5,	6,	2,	4,	30.9771,
5,	6,	2,	3,	37.8575,
5,	6,	2,	2,	41.1933,
5,	6,	2,	1,	39.1733,
5,	6,	2,	0,	39.7539,
5,	6,	1,	4,	34.3563,
5,	6,	1,	3,	42.1185,
5,	6,	1,	2,	45.8597,
5,	6,	1,	1,	43.57,
5,	6,	1,	0,	44.2154,
5,	6,	0,	4,	43.8543,
5,	6,	0,	3,	54.0994,
5,	6,	0,	2,	58.9816,
5,	6,	0,	1,	55.9319,
5,	6,	0,	0,	56.7596,
5,	7,	4,	3,	24.6756,
5,	7,	4,	2,	28.3646,
5,	7,	4,	1,	24.8194,
5,	7,	4,	0,	25.5238,
5,	7,	3,	3,	23.5943,
5,	7,	3,	2,	27.1082,
5,	7,	3,	1,	23.7319,
5,	7,	3,	0,	24.4027,
5,	7,	2,	3,	27.3736,
5,	7,	2,	2,	31.5325,
5,	7,	2,	1,	27.5284,
5,	7,	2,	0,	28.3211,
5,	7,	1,	3,	30.3075,
5,	7,	1,	2,	34.9637,
5,	7,	1,	1,	30.4761,
5,	7,	1,	0,	31.3631,
5,	7,	0,	3,	38.552,
5,	7,	0,	2,	44.6077,
5,	7,	0,	1,	38.7595,
5,	7,	0,	0,	39.9117,
5,	8,	4,	2,	18.3499,
5,	8,	4,	1,	18.8273,
5,	8,	4,	0,	19.5615,
5,	8,	3,	2,	17.5718,
5,	8,	3,	1,	18.0256,
5,	8,	3,	0,	18.7238,
5,	8,	2,	2,	20.2311,
5,	8,	2,	1,	20.7736,
5,	8,	2,	0,	21.607,
5,	8,	1,	2,	22.3016,
5,	8,	1,	1,	22.9122,
5,	8,	1,	0,	23.8497,
5,	8,	0,	2,	28.1164,
5,	8,	0,	1,	28.919,
5,	8,	0,	0,	30.1493,
5,	9,	4,	1,	13.5463,
5,	9,	4,	0,	15.1789,
5,	9,	3,	1,	12.991,
5,	9,	3,	0,	14.5442,
5,	9,	2,	1,	14.8515,
5,	9,	2,	0,	16.7019,
5,	9,	1,	1,	16.3038,
5,	9,	1,	0,	18.3828,
5,	9,	0,	1,	20.3801,
5,	9,	0,	0,	23.1028,
5,	10,	4,	0,	11.9812,
5,	10,	3,	0,	11.491,
5,	10,	2,	0,	13.1417,
5,	10,	1,	0,	14.4292,
5,	10,	0,	0,	18.0433,
5,	15,	4,	3,	79.3176,
5,	15,	4,	2,	94.0521,
5,	15,	4,	1,	78.3857,
5,	15,	4,	0,	80.1693,
5,	15,	3,	3,	75.6345,
5,	15,	3,	2,	89.6313,
5,	15,	3,	1,	74.7472,
5,	15,	3,	0,	76.4433,
5,	15,	2,	3,	89.0831,
5,	15,	2,	2,	106.0054,
5,	15,	2,	1,	88.0314,
5,	15,	2,	0,	90.0631,
5,	15,	1,	3,	99.4643,
5,	15,	1,	2,	118.6201,
5,	15,	1,	1,	98.2857,
5,	15,	1,	0,	100.5748,
5,	15,	0,	3,	128.6652,
5,	15,	0,	2,	154.1105,
5,	15,	0,	1,	127.1296,
5,	15,	0,	0,	130.1436,
5,	16,	4,	2,	74.5165,
5,	16,	4,	1,	71.0924,
5,	16,	4,	0,	68.0049,
5,	16,	3,	2,	71.0934,
5,	16,	3,	1,	67.8356,
5,	16,	3,	0,	64.8988,
5,	16,	2,	2,	83.4671,
5,	16,	2,	1,	79.5813,
5,	16,	2,	0,	76.0713,
5,	16,	1,	2,	93.0311,
5,	16,	1,	1,	88.6627,
5,	16,	1,	0,	84.7125,
5,	16,	0,	2,	119.9288,
5,	16,	0,	1,	114.2016,
5,	16,	0,	0,	109.0125,
5,	17,	4,	1,	62.0011,
5,	17,	4,	0,	57.9021,
5,	17,	3,	1,	59.2049,
5,	17,	3,	0,	55.3028,
5,	17,	2,	1,	69.1518,
5,	17,	2,	0,	64.5162,
5,	17,	1,	1,	76.8564,
5,	17,	1,	0,	71.6558,
5,	17,	0,	1,	98.5175,
5,	17,	0,	0,	91.7271,
5,	18,	4,	0,	49.5774,
5,	18,	3,	0,	47.3889,
5,	18,	2,	0,	55.0478,
5,	18,	1,	0,	60.9928,
5,	18,	0,	0,	77.7004,
5,	35,	4,	1,	86.1092,
5,	35,	4,	0,	79.5494,
5,	35,	3,	1,	82.1821,
5,	35,	3,	0,	75.9392,
5,	35,	2,	1,	96.2974,
5,	35,	2,	0,	88.8643,
5,	35,	1,	1,	107.2157,
5,	35,	1,	0,	98.8671,
5,	35,	0,	1,	137.9185,
5,	35,	0,	0,	126.9932,
6,	6,	4,	4,	23.8856,
6,	6,	4,	3,	28.8016,
6,	6,	4,	2,	31.2455,
6,	6,	4,	1,	29.8333,
6,	6,	4,	0,	30.2736,
6,	6,	3,	3,	35.0224,
6,	6,	3,	2,	38.0624,
6,	6,	3,	1,	36.2493,
6,	6,	3,	0,	36.7837,
6,	6,	2,	2,	41.3845,
6,	6,	2,	1,	39.3923,
6,	6,	2,	0,	39.9744,
6,	6,	1,	1,	37.5254,
6,	6,	1,	0,	38.0804,
6,	6,	0,	0,	38.6444,
6,	7,	4,	3,	21.2786,
6,	7,	4,	2,	24.3496,
6,	7,	4,	1,	21.4046,
6,	7,	4,	0,	21.9902,
6,	7,	3,	3,	25.5404,
6,	7,	3,	2,	29.3414,
6,	7,	3,	1,	25.6857,
6,	7,	3,	0,	26.4095,
6,	7,	2,	3,	27.673,
6,	7,	2,	2,	31.8221,
6,	7,	2,	1,	27.8302,
6,	7,	2,	0,	28.6206,
6,	7,	1,	3,	26.4572,
6,	7,	1,	2,	30.3886,
6,	7,	1,	1,	26.6097,
6,	7,	1,	0,	27.3591,
6,	7,	0,	3,	26.8437,
6,	7,	0,	2,	30.8346,
6,	7,	0,	1,	26.9991,
6,	7,	0,	0,	27.7603,
6,	8,	4,	2,	16.0362,
6,	8,	4,	1,	16.4264,
6,	8,	4,	0,	17.0283,
6,	8,	3,	2,	19.0301,
6,	8,	3,	1,	19.5207,
6,	8,	3,	0,	20.2759,
6,	8,	2,	2,	20.5597,
6,	8,	2,	1,	21.0974,
6,	8,	2,	0,	21.9243,
6,	8,	1,	2,	19.7231,
6,	8,	1,	1,	20.2302,
6,	8,	1,	0,	21.0108,
6,	8,	0,	2,	20.007,
6,	8,	0,	1,	20.5219,
6,	8,	0,	0,	21.3144,
6,	9,	4,	1,	11.9918,
6,	9,	4,	0,	13.3356,
6,	9,	3,	1,	14.0828,
6,	9,	3,	0,	15.7634,
6,	9,	2,	1,	15.1717,
6,	9,	2,	0,	17.0104,
6,	9,	1,	1,	14.5993,
6,	9,	1,	0,	16.3363,
6,	9,	0,	1,	14.8052,
6,	9,	0,	0,	16.5686,
6,	10,	4,	0,	10.6111,
6,	10,	3,	0,	12.4667,
6,	10,	2,	0,	13.4291,
6,	10,	1,	0,	12.9194,
6,	10,	0,	0,	13.1002,
6,	15,	4,	3,	66.6804,
6,	15,	4,	2,	78.6013,
6,	15,	4,	1,	65.9078,
6,	15,	4,	0,	67.3671,
6,	15,	3,	3,	81.8881,
6,	15,	3,	2,	97.1254,
6,	15,	3,	1,	80.9292,
6,	15,	3,	0,	82.769,
6,	15,	2,	3,	89.2013,
6,	15,	2,	2,	105.9214,
6,	15,	2,	1,	88.1536,
6,	15,	2,	0,	90.1686,
6,	15,	1,	3,	84.7049,
6,	15,	1,	2,	100.3938,
6,	15,	1,	1,	83.7129,
6,	15,	1,	0,	85.6116,
6,	15,	0,	3,	85.9634,
6,	15,	0,	2,	101.8755,
6,	15,	0,	1,	84.9562,
6,	15,	0,	0,	86.883,
6,	16,	4,	2,	62.9588,
6,	16,	4,	1,	60.141,
6,	16,	4,	0,	57.607,
6,	16,	3,	2,	76.9439,
6,	16,	3,	1,	73.4146,
6,	16,	3,	0,	70.2311,
6,	16,	2,	2,	83.7301,
6,	16,	2,	1,	79.8688,
6,	16,	2,	0,	76.3841,
6,	16,	1,	2,	79.6254,
6,	16,	1,	1,	75.9798,
6,	16,	1,	0,	72.6928,
6,	16,	0,	2,	80.809,
6,	16,	0,	1,	77.1088,
6,	16,	0,	0,	73.7728,
6,	17,	4,	1,	52.821,
6,	17,	4,	0,	49.429,
6,	17,	3,	1,	64.0534,
6,	17,	3,	0,	59.8305,
6,	17,	2,	1,	69.5828,
6,	17,	2,	0,	64.9677,
6,	17,	1,	1,	66.3265,
6,	17,	1,	0,	61.9611,
6,	17,	0,	1,	67.3106,
6,	17,	0,	0,	62.8792,
6,	18,	4,	0,	42.629,
6,	18,	3,	0,	51.2683,
6,	18,	2,	0,	55.585,
6,	18,	1,	0,	53.1146,
6,	18,	0,	0,	53.8974,
6,	35,	4,	1,	72.9906,
6,	35,	4,	0,	67.5779,
6,	35,	3,	1,	88.9392,
6,	35,	3,	0,	82.1782,
6,	35,	2,	1,	96.7179,
6,	35,	2,	0,	89.3249,
6,	35,	1,	1,	92.0575,
6,	35,	1,	0,	85.0717,
6,	35,	0,	1,	93.424,
6,	35,	0,	0,	86.3333,
7,	7,	3,	3,	19.0289,
7,	7,	3,	2,	21.7161,
7,	7,	3,	1,	19.1404,
7,	7,	3,	0,	19.6515,
7,	7,	2,	2,	24.8343,
7,	7,	2,	1,	21.8429,
7,	7,	2,	0,	22.4366,
7,	7,	1,	1,	19.2533,
7,	7,	1,	0,	19.7676,
7,	7,	0,	0,	20.2979,
7,	8,	3,	2,	14.4583,
7,	8,	3,	1,	14.7957,
7,	8,	3,	0,	15.3174,
7,	8,	2,	2,	16.3999,
7,	8,	2,	1,	16.795,
7,	8,	2,	0,	17.405,
7,	8,	1,	2,	14.5435,
7,	8,	1,	1,	14.8829,
7,	8,	1,	0,	15.4075,
7,	8,	0,	2,	14.9115,
7,	8,	0,	1,	15.262,
7,	8,	0,	0,	15.8036,
7,	9,	3,	1,	10.9012,
7,	9,	3,	0,	12.0692,
7,	9,	2,	1,	12.2933,
7,	9,	2,	0,	13.6561,
7,	9,	1,	1,	10.9648,
7,	9,	1,	0,	12.1393,
7,	9,	0,	1,	11.2272,
7,	9,	0,	0,	12.439,
7,	10,	3,	0,	9.6567,
7,	10,	2,	0,	10.8862,
7,	10,	1,	0,	9.7117,
7,	10,	0,	0,	9.9428,
7,	15,	3,	3,	58.7457,
7,	15,	3,	2,	69.0486,
7,	15,	3,	1,	58.0719,
7,	15,	3,	0,	59.3383,
7,	15,	2,	3,	67.8382,
7,	15,	2,	2,	79.9498,
7,	15,	2,	1,	67.0546,
7,	15,	2,	0,	68.536,
7,	15,	1,	3,	59.073,
7,	15,	1,	2,	69.4142,
7,	15,	1,	1,	58.3952,
7,	15,	1,	0,	59.6678,
7,	15,	0,	3,	60.8063,
7,	15,	0,	2,	71.4864,
7,	15,	0,	1,	60.1073,
7,	15,	0,	0,	61.4207,
7,	16,	3,	2,	55.617,
7,	16,	3,	1,	53.1655,
7,	16,	3,	0,	50.9634,
7,	16,	2,	2,	64.0769,
7,	16,	2,	1,	61.2168,
7,	16,	2,	0,	58.6447,
7,	16,	1,	2,	55.9343,
7,	16,	1,	1,	53.47,
7,	16,	1,	0,	51.2569,
7,	16,	0,	2,	57.5483,
7,	16,	0,	1,	55.0061,
7,	16,	0,	0,	52.7226,
7,	17,	3,	1,	46.8778,
7,	17,	3,	0,	43.9193,
7,	17,	2,	1,	53.8021,
7,	17,	2,	0,	50.3593,
7,	17,	1,	1,	47.153,
7,	17,	1,	0,	44.1781,
7,	17,	0,	1,	48.4747,
7,	17,	0,	0,	45.4072,
7,	18,	3,	0,	38.038,
7,	18,	2,	0,	43.4708,
7,	18,	1,	0,	38.2648,
7,	18,	0,	0,	39.3008,
7,	35,	3,	1,	64.6025,
7,	35,	3,	0,	59.8868,
7,	35,	2,	1,	74.3159,
7,	35,	2,	0,	68.8213,
7,	35,	1,	1,	64.9743,
7,	35,	1,	0,	60.2336,
7,	35,	0,	1,	66.8274,
7,	35,	0,	0,	61.9381,
8,	8,	2,	2,	11.1827,
8,	8,	2,	1,	11.4197,
8,	8,	2,	0,	11.788,
8,	8,	1,	1,	11.6647,
8,	8,	1,	0,	12.0451,
8,	8,	0,	0,	12.444,
8,	9,	2,	1,	8.5761,
8,	9,	2,	0,	9.4059,
8,	9,	1,	1,	8.7411,
8,	9,	1,	0,	9.5974,
8,	9,	0,	1,	8.9987,
8,	9,	0,	0,	9.8956,
8,	10,	2,	0,	7.609,
8,	10,	1,	0,	7.7546,
8,	10,	0,	0,	7.9822,
8,	15,	2,	3,	43.099,
8,	15,	2,	2,	50.2539,
8,	15,	2,	1,	42.6154,
8,	15,	2,	0,	43.5086,
8,	15,	1,	3,	44.2959,
8,	15,	1,	2,	51.7021,
8,	15,	1,	1,	43.7974,
8,	15,	1,	0,	44.7202,
8,	15,	0,	3,	46.1322,
8,	15,	0,	2,	53.9215,
8,	15,	0,	1,	45.6111,
8,	15,	0,	0,	46.5788,
8,	16,	2,	2,	41.0843,
8,	16,	2,	1,	39.3414,
8,	16,	2,	0,	37.7816,
8,	16,	1,	2,	42.1885,
8,	16,	1,	1,	40.39,
8,	16,	1,	0,	38.7797,
8,	16,	0,	2,	43.8847,
8,	16,	0,	1,	42.0013,
8,	16,	0,	0,	40.3138,
8,	17,	2,	1,	35.022,
8,	17,	2,	0,	32.9025,
8,	17,	1,	1,	35.9125,
8,	17,	1,	0,	33.7275,
8,	17,	0,	1,	37.2835,
8,	17,	0,	0,	34.9987,
8,	18,	2,	0,	28.7757,
8,	18,	1,	0,	29.462,
8,	18,	0,	0,	30.5219,
8,	35,	2,	1,	47.9416,
8,	35,	2,	0,	44.5761,
8,	35,	1,	1,	49.2023,
8,	35,	1,	0,	45.7311,
8,	35,	0,	1,	51.1406,
8,	35,	0,	0,	47.5081,
9,	9,	1,	1,	6.6851,
9,	9,	1,	0,	7.2698,
9,	9,	0,	0,	7.9435,
9,	10,	1,	0,	5.9454,
9,	10,	0,	0,	6.4619,
9,	15,	1,	3,	31.4132,
9,	15,	1,	2,	36.3589,
9,	15,	1,	1,	31.069,
9,	15,	1,	0,	31.695,
9,	15,	0,	3,	35.482,
9,	15,	0,	2,	41.2684,
9,	15,	0,	1,	35.088,
9,	15,	0,	0,	35.8129,
9,	16,	1,	2,	30.1398,
9,	16,	1,	1,	28.9095,
9,	16,	1,	0,	27.8121,
9,	16,	0,	2,	33.9048,
9,	16,	0,	1,	32.4877,
9,	16,	0,	0,	31.2207,
9,	17,	1,	1,	25.9699,
9,	17,	1,	0,	24.4631,
9,	17,	0,	1,	29.0227,
9,	17,	0,	0,	27.2958,
9,	18,	1,	0,	21.5952,
9,	18,	0,	0,	23.9644,
9,	35,	1,	1,	35.3301,
9,	35,	1,	0,	32.9451,
9,	35,	0,	1,	39.6376,
9,	35,	0,	0,	36.8974,
10,	10,	0,	0,	5.2971,
10,	15,	0,	3,	27.7982,
10,	15,	0,	2,	32.1989,
10,	15,	0,	1,	27.4945,
10,	15,	0,	0,	28.0493,
10,	16,	0,	2,	26.664,
10,	16,	0,	1,	25.5754,
10,	16,	0,	0,	24.6038,
10,	17,	0,	1,	22.9724,
10,	17,	0,	0,	21.641,
10,	18,	0,	0,	19.1104,
10,	35,	0,	1,	31.2612,
10,	35,	0,	0,	29.1518,
15,	15,	3,	3,	193.8204,
15,	15,	3,	2,	231.3861,
15,	15,	3,	1,	191.5166,
15,	15,	3,	0,	195.9997,
15,	15,	2,	2,	277.4987,
15,	15,	2,	1,	228.6213,
15,	15,	2,	0,	234.0661,
15,	15,	1,	1,	189.2407,
15,	15,	1,	0,	193.669,
15,	15,	0,	0,	198.2084,
15,	16,	3,	2,	181.109,
15,	16,	3,	1,	172.5613,
15,	16,	3,	0,	164.8289,
15,	16,	2,	2,	215.4786,
15,	16,	2,	1,	205.1453,
15,	16,	2,	0,	195.777,
15,	16,	1,	2,	178.9685,
15,	16,	1,	1,	170.5251,
15,	16,	1,	0,	162.8872,
15,	16,	0,	2,	183.0993,
15,	16,	0,	1,	174.4474,
15,	16,	0,	0,	166.6194,
15,	17,	3,	1,	149.3706,
15,	17,	3,	0,	139.2041,
15,	17,	2,	1,	176.7527,
15,	17,	2,	0,	164.5158,
15,	17,	1,	1,	147.6237,
15,	17,	1,	0,	137.5806,
15,	17,	0,	1,	150.9515,
15,	17,	0,	0,	140.6643,
15,	18,	3,	0,	118.3064,
15,	18,	2,	0,	139.1963,
15,	18,	1,	0,	116.9411,
15,	18,	0,	0,	119.5084,
15,	35,	3,	1,	208.5868,
15,	35,	3,	0,	192.2588,
15,	35,	2,	1,	247.6696,
15,	35,	2,	0,	227.9675,
15,	35,	1,	1,	206.1321,
15,	35,	1,	0,	190.0027,
15,	35,	0,	1,	210.8477,
15,	35,	0,	0,	194.3229,
16,	16,	2,	2,	169.6896,
16,	16,	2,	1,	161.7869,
16,	16,	2,	0,	154.6493,
16,	16,	1,	1,	154.2771,
16,	16,	1,	0,	147.4969,
16,	16,	0,	0,	141.0422,
16,	17,	2,	1,	140.5728,
16,	17,	2,	0,	131.1431,
16,	17,	1,	1,	134.1713,
16,	17,	1,	0,	125.2037,
16,	17,	0,	1,	128.4044,
16,	17,	0,	0,	119.8563,
16,	18,	2,	0,	111.8748,
16,	18,	1,	0,	106.908,
16,	18,	0,	0,	102.4455,
16,	35,	2,	1,	195.7624,
16,	35,	2,	0,	180.6446,
16,	35,	1,	1,	186.7232,
16,	35,	1,	0,	172.3521,
16,	35,	0,	1,	178.5666,
16,	35,	0,	0,	164.8741,
16,	17,	1,	1,	117.2993,
16,	17,	1,	0,	109.6211,
16,	17,	0,	0,	102.4888,
16,	18,	1,	0,	94.095,
16,	18,	0,	0,	88.1052,
16,	35,	1,	1,	162.6241,
16,	35,	1,	0,	150.349,
16,	35,	0,	1,	151.8182,
16,	35,	0,	0,	140.4228,
18,	18,	0,	0,	76.145,
18,	35,	0,	1,	129.8276,
18,	35,	0,	0,	120.2787,
35,	35,	1,	1,	226.0879,
35,	35,	1,	0,	208.7808,
35,	35,	0,	0,	192.8938
};
  //  cout <<  C6_Matrix[4] << endl;
    int Ntot = GetTotalNumberOfAtoms();
    double *p;
    p = Coord_Num::coord_num().CalculateCN2();


} */






















