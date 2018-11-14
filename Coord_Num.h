#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
//#include"Dispersion.h"

using namespace std;

class Coord_Num {
    private:
        Coord_Num();
        ~Coord_Num();

    double *CutoffRadii, *r0AB, *Coordination_Number, *Coordination_Number_step, *Coordination_Number_Deriv,*Coordination_Number_FD, *UCHF, *CKS, *GaussA_C6_CKS,*C6_CKS_deriv,*C6_CKS_FD, *d_E6_CKS, *modR_AB, *GaussA_C6_UCHF, *C6_UCHF_deriv, *C6_UCHF_FD, *d_E6_UCHF, *C8CKS, *C8CKS_deriv, *d_E8_CKS, *C8UCHF,  *C8UCHF_deriv, *d_E8_UCHF, *E6_CKS_coord_derivs, *E6_UCHF_coord_derivs, *E8_CKS_coord_derivs, *E8_UCHF_coord_derivs; //*RcovAB; */
    double **RcovAB;
    int *Ref_Compounds;
    //double *CovalentRadii;

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
       // int GetCovalentRadii();
        // formerly this was GetCovalentRadii
        int GetCutoffRadii();
        int GetCovalentRadii();
        int GetCovalentRadii2();
        int CalculateCN();
        //int CalculateCN2();
        //int CalculateCN2();
        double *CalculateCN2();
        double *CalculateCN2_step();
        double *CalculateCN2_derivative();
        double *CalculateCN2_derivative_FD();
        int *GetNumberReferenceCompounds();
//        double GetC6();
        double GetCKSC6();
        double GetUCHFC6();
        double GetCKSEnergy();
        double GetCKSC8Energy();
        double GetUCHFEnergy();
        double GetUCHFC8Energy();
        double MP2DDispersionCorrection();
        double GetGradient();
        
        
        
};
