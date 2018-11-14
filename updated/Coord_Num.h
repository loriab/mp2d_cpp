#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>

using namespace std;

class Coord_Num {
    private:
        Coord_Num();
        ~Coord_Num();

    double *CutoffRadii, *r0AB, *Coordination_Number, *Coordination_Number_step, *Coordination_Number_Deriv,*Coordination_Number_FD, *UCHF, *CKS, *GaussA_C6_CKS,*C6_CKS_deriv,*C6_CKS_FD, *d_E6_CKS, *modR_AB, *GaussA_C6_UCHF, *C6_UCHF_deriv, *C6_UCHF_FD, *d_E6_UCHF, *C8CKS, *C8CKS_deriv, *d_E8_CKS, *C8UCHF,  *C8UCHF_deriv, *d_E8_UCHF, *E6_CKS_coord_derivs, *E6_UCHF_coord_derivs, *E8_CKS_coord_derivs, *E8_UCHF_coord_derivs;
    double **RcovAB;
    int *Ref_Compounds;
    int C6_counter;
    int UCHF_counter;
    double a_one, a_two, rcut, width, s_8;
    vector<vector<double> > matrix_Grimme;
    vector<vector<double> > matrix_UCHF;

    public:
        // Instantiate the single instance of the CN class
        static Coord_Num& coord_num() {
            static Coord_Num thecoord_num;
            return thecoord_num;
        }

        double L_ij;
        double L_ij_step;
        double L_ij_deriv;
        double deriv_prefactor;
        double deriv_postfactor;

        string Path;
        string GrimmePath;
        string UCHFPath;
        
        

        void Initialize(ifstream& infile);
        void UseFunction();
        // Can I grab some of the parameters from the input file?
        void GetParameters(ifstream& infile);
        void GetCutoffRadii();
        int GetCovalentRadii();
        void GetCovalentRadii2();
        int CalculateCN();
        double *CalculateCN2();
        double *CalculateCN2_step();
        double *CalculateCN2_derivative();
        double *CalculateCN2_derivative_FD();
        void *GetNumberReferenceCompounds();
        void GetC6Coefficients();
        void GetCKSC6();
        void GetUCHFC6();
        double GetCKSEnergy();
        double GetCKSC8Energy();
        double GetUCHFEnergy();
        double GetUCHFC8Energy();
        void MP2DDispersionCorrection();
        void GetGradient();
        
        void Rewind(ifstream& infile) {
        infile.clear();
        infile.seekg(0,ios::beg);
        }
        
        
};
