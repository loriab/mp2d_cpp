#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"Dispersion.h"
#include"Coord_Num.h"
#include<iterator>
#include<unistd.h>



//Constructor
Coord_Num::Coord_Num() : /*CovalentRadii(NULL),*/ CutoffRadii(NULL),r0AB(NULL), RcovAB(NULL), Coordination_Number(NULL), Coordination_Number_step(NULL), Coordination_Number_Deriv(NULL),Coordination_Number_FD(NULL), Ref_Compounds(NULL), UCHF(NULL), CKS(NULL), GaussA_C6_CKS(NULL),C6_CKS_deriv(NULL), C6_CKS_FD(NULL), d_E6_CKS(NULL), E6_CKS_coord_derivs(NULL), modR_AB(NULL), GaussA_C6_UCHF(NULL), C6_UCHF_deriv(NULL), C6_UCHF_FD(NULL), d_E6_UCHF(NULL),  E6_UCHF_coord_derivs(NULL), C8CKS(NULL), C8CKS_deriv(NULL), E8_CKS_coord_derivs(NULL), d_E8_CKS(NULL), C8UCHF(NULL), C8UCHF_deriv(NULL),  E8_UCHF_coord_derivs(NULL), d_E8_UCHF(NULL)
{
}

// Need to initialize functions
void Coord_Num::Initialize(ifstream &infile) {
//    UseFunction();
//    GetCovalentRadii();
    // formerly this was GetCovalentRadii
    GetCutoffRadii();
//    GetCovalentRadii();
    GetCovalentRadii2();
//    CalculateCN();
    CalculateCN2();
    CalculateCN2_step();
    CalculateCN2_derivative();
    CalculateCN2_derivative_FD();
    GetNumberReferenceCompounds();
//    GetC6();
    GetCKSC6();
    GetUCHFC6();
    GetCKSEnergy();
    GetCKSC8Energy();
    GetUCHFEnergy();
    GetUCHFC8Energy();
    MP2DDispersionCorrection();
    GetGradient();
    
  
}

//Destructor
Coord_Num::~Coord_Num() {

    delete [] CutoffRadii;
    delete [] r0AB;
    delete [] RcovAB;
    delete [] Coordination_Number;
    delete [] Coordination_Number_step;
    delete [] Coordination_Number_Deriv;
    delete [] Coordination_Number_FD;
    delete [] Ref_Compounds;


    delete [] CKS;
    delete [] UCHF;

    delete [] GaussA_C6_CKS;
    delete [] C6_CKS_deriv;
    delete [] C6_CKS_FD;
    delete [] GaussA_C6_UCHF;
    delete [] C6_UCHF_deriv;
    delete [] C6_UCHF_FD;
    delete [] C8CKS;
    delete [] C8CKS_deriv;
    delete [] C8UCHF;
    delete [] C8UCHF_deriv;
 //   delete [] E8_CKS_coord_derivs;
    delete [] E8_UCHF_coord_derivs;
    
 
}

// This is just a test function I used to learn how to call functions from other classes.
/*void Coord_Num::UseFunction() {
    int Ntot;
   // Ntot = Dispersion::GetTotalNumberOfAtoms();
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();
    cout << Ntot << "from the Coordination Number Class" << endl;
}*/

//int Coord_Num::GetCovalentRadii() {
int Coord_Num::GetCutoffRadii() {

    // an array that contains the r0AB values for pairs of atoms.
    // For example, CovalentRadii[1][1] is the r0AB for HH.
    // CovalentRadii[1][5] is the r0AB for HB.
//    double CovalentRadii[36][36] = {
    double CutoffRadii[36][36] = {
        0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	2.1823,	0.0000,	0.0000,	0.0000,	2.5141,	2.4492,	2.3667,	2.1768,	2.0646,	1.9892,	0.0000,	0.0000,	0.0000,	0.0000,	2.8304,	2.6190,	2.4757,	2.3725,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.6026,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	2.5141,	0.0000,	0.0000,	0.0000,	3.2160,	2.9531,	2.7776,	2.6482,	2.6233,	2.4994,	0.0000,	0.0000,	0.0000,	0.0000,	3.1866,	3.0651,	2.9768,	2.9093,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.1029,
	0.0000,	2.4492,	0.0000,	0.0000,	0.0000,	2.9531,	2.9103,	2.7063,	2.5697,	2.4770,	2.4091,	0.0000,	0.0000,	0.0000,	0.0000,	3.1245,	2.9879,	2.8848,	2.8040,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.0108,
	0.0000,	2.3667,	0.0000,	0.0000,	0.0000,	2.7776,	2.7063,	2.6225,	2.4846,	2.3885,	2.3176,	0.0000,	0.0000,	0.0000,	0.0000,	3.0465,	2.9054,	2.7952,	2.7071,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.9227,
	0.0000,	2.1768,	0.0000,	0.0000,	0.0000,	2.6482,	2.5697,	2.4846,	2.4817,	2.3511,	2.2571,	0.0000,	0.0000,	0.0000,	0.0000,	2.8727,	2.8805,	2.7457,	2.6386,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.8694,
	0.0000,	2.0646,	0.0000,	0.0000,	0.0000,	2.6233,	2.4770,	2.3885,	2.3511,	2.2996,	2.1946,	0.0000,	0.0000,	0.0000,	0.0000,	2.7664,	2.7330,	2.6881,	2.5720,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.8109,
	0.0000,	1.9892,	0.0000,	0.0000,	0.0000,	2.4994,	2.4091,	2.3176,	2.2571,	2.1946,	2.1374,	0.0000,	0.0000,	0.0000,	0.0000,	2.6926,	2.6331,	2.5728,	2.5139,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.6929,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	2.8304,	0.0000,	0.0000,	0.0000,	3.1866,	3.1245,	3.0465,	2.8727,	2.7664,	2.6926,	0.0000,	0.0000,	0.0000,	0.0000,	3.5017,	3.3180,	3.1916,	3.0982,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.3218,
	0.0000,	2.6190,	0.0000,	0.0000,	0.0000,	3.0651,	2.9879,	2.9054,	2.8805,	2.7330,	2.6331,	0.0000,	0.0000,	0.0000,	0.0000,	3.3180,	3.3107,	3.1523,	3.0352,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.2815,
	0.0000,	2.4757,	0.0000,	0.0000,	0.0000,	2.9768,	2.8848,	2.7952,	2.7457,	2.6881,	2.5728,	0.0000,	0.0000,	0.0000,	0.0000,	3.1916,	3.1523,	3.1046,	2.9730,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.2346,
	0.0000,	2.3725,	0.0000,	0.0000,	0.0000,	2.9093,	2.8040,	2.7071,	2.6386,	2.5720,	2.5139,	0.0000,	0.0000,	0.0000,	0.0000,	3.0982,	3.0352,	2.9730,	2.9148,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.0994,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
	0.0000,	2.6026,	0.0000,	0.0000,	0.0000,	3.1029,	3.0108,	2.9227,	2.8694,	2.8109,	2.6929,	0.0000,	0.0000,	0.0000,	0.0000,	3.3218,	3.2815,	3.2346,	3.0994,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.3662,

    };
    // Just a test print statement to see if the above matrix gives me the r0AB values I want
    /*    cout << CovalentRadii[1][1] << endl;
    cout << CovalentRadii[1][9] << endl;
    cout << CovalentRadii[9][1] << endl; */
    int Ntot;
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int r0ABtot =0;
    for (int i = 0; i<Ntot; i++) {
        r0ABtot += i;
    }
//    cout << Ntot << endl;
    
    // Here is where things get schwifty.
    // I define a pointer *p. I set p equal to the argument of GetAtomicNumber which is a pointer to the array atomic_number.
    // I cycle from 0 to Ntot -1.Then, *(p + i) is the ith element of atomic_number[i].

    int *p;
    p = Dispersion::dispersion().GetAtomicNumber();
    int index=-1;
//    double r0AB[r0ABtot];
    r0AB = new double[r0ABtot];
        for (int i=0;i<Ntot;i++) {
            for (int j = i+1; j<Ntot; j++) {
    //        cout << *(p + i) << endl;
               // cout << CovalentRadii[*(p+i)][*(p+j)] << endl;
                index++;
//                cout << index << endl; 
//                r0AB[index] = CovalentRadii[*(p+i)][*(p+j)];
                r0AB[index] = CutoffRadii[*(p+i)][*(p+j)];
            }
        }
    // A simple print statement loop to check if I'm printing out the correct
    // r0AB for each atom pair AB
/*        for (int i=0;i<r0ABtot;i++) {
            cout << r0AB[i] << endl;
        }*/
}


int Coord_Num::GetCovalentRadii2() {
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

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int *p;
    p = Dispersion::dispersion().GetAtomicNumber();

    int RcovTot=0;
    for (int i=0;i<Ntot;i++) {
        RcovTot += i;
    }
//    cout << RcovTot << endl;

    double Rcov[Ntot];
    for (int i=0; i<Ntot; i++ ) {
//        cout << CovalentRadii[atomic_number[i]] << endl;
        Rcov[i] = CovalentRadii[*(p+i)];
//        cout << CovalentRadii[*(p+i)] << endl;
//        cout << Rcov[i] << endl;
    }

//    double RcovAB[RcovTot];
//    RcovAB = new double[RcovTot];
  //  double RcovAB[Ntot][Ntot];
//    int index = -1;
    RcovAB = new double*[Ntot];
    for(int i=0; i < Ntot; i++) {
        RcovAB[i] = new double[Ntot];
    }
    for (int i=0; i<Ntot; i++) {
        for (int j= 0; j<Ntot; j++) {
           // index++;
            if (i==j) {
                RcovAB[i][j] = 0.0;
            }
            else {
                RcovAB[i][j] = Rcov[i] + Rcov[j];
             //   cout << RcovAB[index] << endl;
            }
//        cout << "the covalent radii is " << RcovAB[i][j] << endl;
        }
    //cout << "the covalent radii is " << RcovAB[i][j] << endl;
    }
//    double test = std::begin(RcovAB);
//    cout << test << endl;
}


double *Coord_Num::CalculateCN2() {

    int Ntot;
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();
    double function;
    double CNum=0;
//    double Coordination_Number[Ntot-1];
 //   Coordination_Number = new double[Ntot-1];
    Coordination_Number = new double[Ntot];
    double xprime;
    for (int i =0; i<Ntot; i++) {
        for (int j=0; j<Ntot; j++) {
            if (i==j) {
                function =0;
            }
//            cout << "the interatomic distance is " << * ( *(p+i)+j) << endl;
 //           cout << "the sum of the covalent radii is " << RcovAB[i][j] << endl;

//////Here is the new way to calculate CN
            else {
                double rAB_t = (* ( *(p+i)+j));
             //   if (* ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]) {
                if ( rAB_t <= 0.95*RcovAB[i][j]) {
//                    * ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]; 
                    function = 1.0;
                }
               // else if (* ( *(p+i)+j)) >= 1.75*RcovAB[i][j]) {
                else if (rAB_t >= 1.75*RcovAB[i][j]) {
//                     * ( *(p+i)+j)) >= 1.75*RcovAB[i][j];
                     function = 0.0;
                }
                else {
                    double top = (* ( *(p+i)+j) - 0.95*RcovAB[i][j]);
                    xprime = top/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
                    function = 1.0 - ((-20*pow(xprime, 7.0)) + (70*pow(xprime, 6.0)) - (84*pow(xprime,5.0)) + (35*pow(xprime,4.0)));
                }
                CNum += function;
                if ((j==Ntot-1) || (i==Ntot-1 && j==Ntot-2)) {
                    Coordination_Number[i] = CNum;
             //       cout << "the coordination number is " << Coordination_Number[i] << endl;
                    CNum =0;
                }
            }
        }
    }
    //return Coordination_Number;
//    cout << Ntot << endl;
//  for (int i=0;i<Ntot;i++) {
        // Originally I printed out the actual coordiantion number 
//        cout << Coordination_Number[i] << endl;
        // For now I'm rounding the coord numbers to test out getting C6
//        cout << round(Coordination_Number[i]) << endl;
//    }
    return Coordination_Number;
}

// This is only here so that I can calculate the derivative of C6 by finite difference

double *Coord_Num::CalculateCN2_step() {

    int Ntot;
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();
    double function;
    double CNum=0;
//    double Coordination_Number[Ntot-1];
 //   Coordination_Number = new double[Ntot-1];
    Coordination_Number_step = new double[Ntot];
    double xprime;
    for (int i =0; i<Ntot; i++) {
        for (int j=0; j<Ntot; j++) {
            if (i==j) {
                function =0;
            }
//            cout << "the interatomic distance is " << * ( *(p+i)+j) << endl;
 //           cout << "the sum of the covalent radii is " << RcovAB[i][j] << endl;

//////Here is the new way to calculate CN
            else {
                double rAB_t = (* ( *(p+i)+j))+0.0000001;
             //   if (* ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]) {
                if ( rAB_t <= 0.95*RcovAB[i][j]) {
//                    * ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]; 
                    function = 1.0;
                }
               // else if (* ( *(p+i)+j)) >= 1.75*RcovAB[i][j]) {
                else if (rAB_t >= 1.75*RcovAB[i][j]) {
//                     * ( *(p+i)+j)) >= 1.75*RcovAB[i][j];
                     function = 0.0;
                }
                else {
                    double top = ((* ( *(p+i)+j)+ 0.0000001) - 0.95*RcovAB[i][j]);
                    xprime = top/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
                    function = 1.0 - ((-20*pow(xprime, 7.0)) + (70*pow(xprime, 6.0)) - (84*pow(xprime,5.0)) + (35*pow(xprime,4.0)));
                }
                CNum += function;
                if ((j==Ntot-1) || (i==Ntot-1 && j==Ntot-2)) {
                    Coordination_Number_step[i] = CNum;
          //          cout << "the coordination number step is " << Coordination_Number_step[i] << endl;
                    CNum =0;
                }
            }
        }
    }

    return Coordination_Number_step;
}

double *Coord_Num::CalculateCN2_derivative() {
    int Ntot;
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();
    double function;
    double CNum_prime=0;

    Coordination_Number_Deriv = new double[Ntot];
    double xprime;
    for (int i =0; i<Ntot; i++) {
        for (int j=0; j<Ntot; j++) {
            if (i==j) {
                function =0;
            }
            else {
                double rAB_t = (* ( *(p+i)+j));
             //   if (* ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]) {
                if ( rAB_t <= 0.95*RcovAB[i][j]) {
//                    * ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]; 
                    function = 0.0;
                }
               // else if (* ( *(p+i)+j)) >= 1.75*RcovAB[i][j]) {
                else if (rAB_t >= 1.75*RcovAB[i][j]) {
//                     * ( *(p+i)+j)) >= 1.75*RcovAB[i][j];
                     function = 0.0;
                }
                else {
                    double top = (* ( *(p+i)+j) - 0.95*RcovAB[i][j]);
                    xprime = top/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
                    function = (-140*(pow(top,3)*pow(1.75*RcovAB[i][j]-(*(*(p+i)+j)),3)))/(pow(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j],7));
                }
                CNum_prime += function;
                if ((j==Ntot-1) || (i==Ntot-1 && j==Ntot-2)) {
                    Coordination_Number_Deriv[i] = CNum_prime;
              //      cout << "the derivative of CN is " << Coordination_Number_Deriv[i] << endl;
                    CNum_prime =0;
                }
            }
        }
    }
    return Coordination_Number_Deriv;

}

double *Coord_Num::CalculateCN2_derivative_FD() {
    int Ntot;
    Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();
    double function;
    double CNum_FD=0;

    Coordination_Number_FD = new double[Ntot];
    double xprime;
    double xprime_FD;
    for (int i =0; i<Ntot; i++) {
        for (int j=0; j<Ntot; j++) {
            if (i==j) {
                function =0;
            }
            else {
                double rAB_t = (* ( *(p+i)+j));
             //   if (* ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]) {
                if ( rAB_t <= 0.95*RcovAB[i][j]) {
//                    * ( *(p+i)+j))  <= 0.95 *RcovAB[i][j]; 
                    function = 0.0;
                }
               // else if (* ( *(p+i)+j)) >= 1.75*RcovAB[i][j]) {
                else if (rAB_t >= 1.75*RcovAB[i][j]) {
//                     * ( *(p+i)+j)) >= 1.75*RcovAB[i][j];
                     function = 0.0;
                }
                else {
                    double top_FD = (* ( *(p+i)+j)+0.00001 - 0.95*RcovAB[i][j]);
                    double top = (* ( *(p+i)+j) - 0.95*RcovAB[i][j]);
                    xprime_FD = top_FD/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
                    xprime = top/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
                    function = ((1.0 - ((-20*pow(xprime_FD, 7.0)) + (70*pow(xprime_FD, 6.0)) - (84*pow(xprime_FD,5.0)) + (35*pow(xprime_FD,4.0)))) - (1.0 - ((-20*pow(xprime, 7.0)) + (70*pow(xprime, 6.0)) - (84*pow(xprime,5.0)) + (35*pow(xprime,4.0)))))/(0.00001);
                    
                }
                CNum_FD += function;
                if ((j==Ntot-1) || (i==Ntot-1 && j==Ntot-2)) {
                    Coordination_Number_FD[i] = CNum_FD;
             //       cout << "the derivative of CN by finite difference is " << Coordination_Number_FD[i] << endl;
                    CNum_FD =0;
                }
            }
        }
    }
    return Coordination_Number_FD;

}


int *Coord_Num::GetNumberReferenceCompounds() {
    // I need to know the number of reference compounds per atom pair later on for the derivative of C6 with
    // respect to R
    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();
    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    int *p;
    p = Dispersion::dispersion().GetAtomicNumber();

  //  int MAX_PATH_LEN = 500;
  //  char cwd[MAX_PATH_LEN];
  //  getcwd(cwd, MAX_PATH_LEN);
  //  GrimmePath = cwd;
  //  GrimmePath += "/Grimme_C6.txt";
  //  cout << GrimmePath << endl;
   

    ifstream inFile2;
    int index = -1;
    int counter = 0;
    Ref_Compounds = new int[combinations];

    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++ ) {
            double a = *(p+i);
            double b = *(p+j);
            index++;
         //   inFile2.open("Grimme_C6.txt");
         //   inFile2.open(GrimmePath.c_str());
          //  inFile2.open("home/chandemonium/MP2D/Dispersion/Grimme_C6.txt");
          //  inFile2.open("/home/chandemonium/MP2D/Dispersion/Grimme_C6.txt");
            inFile2.open("/home/chandemonium/MP2D_Program/Dispersion/GrimmeC6.txt");
            string line2;
            if (!inFile2) {
                cerr << "unable to open file Grimme_C6.txt" << endl;
                exit (1);
            }
            while (getline(inFile2, line2)) {
                istringstream ss(line2);
                double var6, var7, var8, var9, var10;
                ss >> var6 >> var7 >> var8 >> var9 >> var10;
                if ((var6 == a && var7 == b) || (var6==b && var7==a)) {
                counter++;
                }
            }
            Ref_Compounds[index] = counter;
    //        cout << Ref_Compounds[index] << endl;
            counter =0;
            inFile2.close();
        }
    }
    for (int i=0; i<combinations;i++) {
 //       cout << "the number of reference compounds is " << Ref_Compounds[i] << endl;
    }
}

double Coord_Num::GetCKSC6() {
  int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();
   // cout << Ntot << endl;

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    int *p;
    p = Dispersion::dispersion().GetAtomicNumber();

    ifstream inFile2;
    

    int index = -1;
 //   double UCHF[combinations];
    CKS = new double[combinations];
//    double GaussA_C6[combinations];
    GaussA_C6_CKS = new double[combinations];
    C6_CKS_deriv = new double[combinations];
    C6_CKS_FD = new double[combinations];
    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++ ) {
            index++;
            // These are used to get the precomputed C6
            double a = *(p+i);
            double b = *(p+j);
            double c = round(Coordination_Number[i]);
            double d = round(Coordination_Number[j]);
            // These are used to get the Gaussian Averaged C6
            double e = Coordination_Number[i];
            double g = Coordination_Number[j];
            double h = Coordination_Number_Deriv[i];
            double k = Coordination_Number_Deriv[j];
        //    double h = Coordination_Number_FD[i];
        //    double k = Coordination_Number_FD[j];

        // for computing the derivative of C6 by finite difference
            double r = Coordination_Number_step[i];
            double t = Coordination_Number_step[j];


//inFile2.open("Grimme_C6.txt");
inFile2.open("/home/chandemonium/MP2D_Program/Dispersion/GrimmeC6.txt");
            string line2;
            if (!inFile2) {
                cerr << "unable to open file Grimme_C6.txt" << endl;
                exit (1);
            }
            // In the while loop I will try to get the guassian averaged C6 along with the precomputed value
            double GA_C6 =0;
            double sum_L_ij=0;
            double sum_C6ref_Lij =0;
            double sum_L_ij_deriv=0;
            double sum_C6ref_Lij_deriv=0;
            double sum_numerator=0;
            double sum_C6_deriv=0;

            // for finite difference 
            double sum_L_ij_step=0;
            double sum_C6ref_Lij_step=0;
            
            
            while (getline(inFile2, line2)) {
                istringstream ss(line2);
                double var6, var7, var8, var9, var10;
                ss >> var6 >> var7 >> var8 >> var9 >> var10;
                // Precomputed CKS C6
                if ((var6 == a && var7 == b && round(var8) == c && round(var9) == d) || (var6==b && var7==a && round(var8)==d && round(var9)==c)) {
                //    cout << var5 << endl;
                CKS[index] = var10;
        //        cout << var10 << endl;
                } 
                // Gaussian Averaged C6
                if ((var6 == a && var7 == b) || (var6==b && var7==a)) {
              //  double L_ij = 
               //     cout << var8 << endl;
               //     cout << var9 << endl;
               //     cout << var10 << endl;
               //     double L_ij = exp(-4*(pow((e - round(var8)), 2) + pow((g - round(var9)), 2)));
                   
                    // unless the coordination numbers are equal, the contribution needs to be counted twice.
                    if (a==b && var8 != var9) {

                     //   L_ij = 2*exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9)))));    
                      //  L_ij = exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9))))) + exp(-4*(((e - round(var9))*(e - round(var9))) + ((g - round(var8))*(g - round(var8)))));
                        // The coordination numbers (hybridization numbers are not rounded to the nearest integer when calculating the Gaussian Averaged Dispersion Coefficients.


                        L_ij = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9)))) + exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))));

                        // for finite difference
                        L_ij_step = exp(-4*(((r - var8)*(r - var8)) + ((t - var9)*(t - var9)))) + exp(-4*(((r - var9)*(r - var9)) + ((t - var8)*(t - var8))));

              //          L_ij_deriv = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))))*log(exp(1))*(-8*((e-var8)*h)+((g-var9)*k)) + exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*log(exp(1))*(-8*((e-var9)*k)+((g-var8)*h));

                     //*   L_ij_deriv = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))))*(-8*(((e-var8)*h)+((g-var9)*k))) + exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*(-8*(((e-var9)*k)+((g-var8)*h)));

                        L_ij_deriv = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))))*(-8*(((e-var8)*h)+((g-var9)*k))) + exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*(-8*(((e-var9)*h)+((g-var8)*k)));
                   

                   //     deriv_prefactor = exp(4*(((e-var8)*(e-var8)) + ((g-var9)*(g-var9)))) + exp(4*(((e-var9)*(e-var9)) + ((g-var8)*(g-var8))));

                //        deriv_postfactor = ((8*(e-var8)*h)+(8*(g-var9)*k)) + ((8*(e-var9)*k)+(8*(g-var8)*h));

                   //     L_ij_deriv = (exp(-4*((e-var8)*(e-var8) + (g-var9)*(g-var9))))*(-8*((e-var8)*h + (g-var9)*k)) + (exp(-4*((e-var9)*(e-var9) + (g-var8)*(g-var8))))*(-8*((e-var9)*k + (g-var8)*h));

              //          V_deriv = ((exp(4*((e-var8)*(e-var8) + (g-var9)*(g-var9))))*(8*((e-var8)*h + (g-var9)*k)) + (exp(4*((e-var9)*(e-var9) + (g-var8)*(g-var8))))*(8*((e-var9)*k + (g-var8)*h)));
                    } 
                    
                    // I have to use to else if statements to change the coordination numbers according to which atom comes first.
                    else if (var6 ==a && var7 ==b) {
                     //   L_ij = exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9)))));
                        L_ij = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))));

                        // for finite difference
                        L_ij_step = exp(-4*(((r - var8)*(r - var8)) + ((t - var9)*(t - var9))));
                   
                //        L_ij_deriv = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))))*log(exp(1))*(-8*((e-var8)*h)+((g-var9)*k));

                        L_ij_deriv = exp(-4*(((e - var8)*(e - var8)) + ((g - var9)*(g - var9))))*(-8*(((e-var8)*h)+((g-var9)*k)));

                //        deriv_prefactor = exp(4*(((e-var8)*(e-var8)) + ((g-var9)*(g-var9))));

                //        deriv_postfactor = ((8*(e-var8)*h)+(8*(g-var9)*k));

             //           L_ij_deriv = (exp(-4*((e-var8)*(e-var8) + (g-var9)*(g-var9))))*(-8*((e-var8)*h + (g-var9)*k)); 

               //         V_deriv = ((exp(4*((e-var8)*(e-var8) + (g-var9)*(g-var9))))*(8*((e-var8)*h + (g-var9)*k)));
                    }
                    else if (var6==b && var7 ==a) {
                     //   L_ij = exp(-4*(((e - round(var9))*(e - round(var9))) + ((g - round(var8))*(g - round(var8)))));
                        L_ij = exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))));

                        // for finite difference
                        L_ij_step = exp(-4*(((r - var9)*(r - var9)) + ((t - var8)*(t - var8))));

         //               L_ij_deriv = exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*log(exp(1))*(-8*((e-var9)*k)+((g-var8)*h));

               //*         L_ij_deriv = exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*(-8*(((e-var9)*k)+((g-var8)*h)));

                          L_ij_deriv = exp(-4*(((e - var9)*(e - var9)) + ((g - var8)*(g - var8))))*(-8*(((e-var9)*h)+((g-var8)*k)));

             //           deriv_prefactor = exp(4*(((e-var9)*(e-var9)) + ((g-var8)*(g-var8))));
 
             //           deriv_postfactor = ((8*(e-var9)*k)+(8*(g-var8)*h));

               //           L_ij_deriv = (exp(-4*((e-var9)*(e-var9) + (g-var8)*(g-var8))))*(-8*((e-var9)*k + (g-var8)*h));

                 //        V_deriv = ((exp(4*((e-var9)*(e-var9) + (g-var8)*(g-var8))))*(8*((e-var9)*k + (g-var8)*h)));
                    }
                    else {}

                    sum_L_ij += L_ij;
                    sum_C6ref_Lij += L_ij*var10;
         //           cout << "sum_L_ij " << sum_L_ij << endl;
                    
                    
                    sum_L_ij_deriv += L_ij_deriv;
                    sum_C6ref_Lij_deriv += L_ij_deriv*var10;
              //      cout << "L_ij derivative" << sum_L_ij_deriv << endl;


                   // for finite difference
                   sum_L_ij_step += L_ij_step;
                   sum_C6ref_Lij_step += L_ij_step*var10;
        //           cout << "sum_L_ij step is " << sum_L_ij_step << endl;

         //           double C6_deriv = (deriv_prefactor*var10*L_ij_deriv + deriv_prefactor*log10(exp(1))*var10*L_ij)/(Ref_Compounds[i]*Ref_Compounds[j]);


            //        double C6_deriv = (deriv_prefactor*var10*L_ij_deriv + deriv_prefactor*log10(exp(1))*var10*L_ij*deriv_postfactor)/(Ref_Compounds[i]*Ref_Compounds[j]);

                //    double C6_deriv = (deriv_prefactor*var10*L_ij_deriv + deriv_prefactor*var10*L_ij*deriv_postfactor)/(Ref_Compounds[i]*Ref_Compounds[j]);

          ////     //     double C6_deriv = (deriv_prefactor*var10*L_ij_deriv)/(Ref_Compounds[i]*Ref_Compounds[j]) + (deriv_prefactor*var10*L_ij*deriv_postfactor)/(Ref_Compounds[i]*Ref_Compounds[j]);
       
                      

       
              //      sum_C6_deriv += C6_deriv;
         //*           sum_C6_deriv += (sum_C6ref_Lij_deriv*sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv)/pow(sum_L_ij,2);
                    

                    double numerator = L_ij*var10;
                    sum_numerator += numerator;

                    // for finite difference
                    


           //         cout << "sum_numerator = " << sum_numerator << endl;
               //     cout << "e = " << e << endl;
                //    cout << "var8 = " << var8 << endl;
                 //   cout << "g = " << g << endl;
                 //   cout << "var9 = " << var9 << endl;
            // Neither equation for GA_C6 is correct as of 4/6/2018
                 //   GA_C6 += ((CKS[index] * L_ij) / L_ij);
                 //   GA_C6 += ((var10 * L_ij) / L_ij);
                 //   cout << "GA_C6 = " << GA_C6 << endl;
                GaussA_C6_CKS[index] = sum_numerator/sum_L_ij;
                
    //            C6_CKS_deriv[index] = sum_C6_deriv;
//This is my working derivative
                   C6_CKS_deriv[index] = (sum_C6ref_Lij_deriv*sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv)/pow(sum_L_ij,2);
 // this is the same as the above
         //       C6_CKS_deriv[index] = (sum_C6ref_Lij_deriv/sum_L_ij) - ((sum_L_ij_deriv*sum_C6ref_Lij)/pow(sum_L_ij,2));

                C6_CKS_FD[index] = ((sum_C6ref_Lij_step/sum_L_ij_step) - (sum_numerator/sum_L_ij))/0.0000001;

                
        //        cout << C6_CKS_deriv[index] << endl;
                }
             //   cout << counter << endl;
             //   cout << "GA_C6 = " << GA_C6 << endl; 
                
            }
        //    cout << " Finite Difference Derivative of C6 CKS " << C6_CKS_FD[index] << endl;
       //     cout << "C6 CKS derivative " << C6_CKS_deriv[index] << endl;
            inFile2.close();
        }
    }
    for (int i=0; i<combinations;i++) {
    //    cout << "CKS = " << CKS[i] << endl;
 //       cout << "CKS = " << GaussA_C6_CKS[i] << endl;
 //           cout << "the derivative of the CKS C6 Coefficient is " << C6_CKS_deriv[i] << endl;
    }  

}


double Coord_Num::GetUCHFC6() {
  int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    int *p;
    p = Dispersion::dispersion().GetAtomicNumber();

  //  int MAX_PATH_LEN = 500;
  //  char cwd[MAX_PATH_LEN];
  //  getcwd(cwd, MAX_PATH_LEN);
  //  UCHFPath = cwd;
  //  UCHFPath += "/UCHF_daug_C6.txt";
  //  cout << UCHFPath << endl;

    ifstream inFile;

    int index = -1;
 //   double UCHF[combinations];
    UCHF = new double[combinations];

    GaussA_C6_UCHF = new double[combinations];

    C6_UCHF_deriv = new double[combinations];

    C6_UCHF_FD = new double[combinations];

    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++ ) {
            index++;
            double a = *(p+i);
            double b = *(p+j);
            double c = round(Coordination_Number[i]);
            double d = round(Coordination_Number[j]);
            double e = Coordination_Number[i];
            double g = Coordination_Number[j];
            // for C6 derivative
            double h = Coordination_Number_Deriv[i];
            double k = Coordination_Number_Deriv[j];

            // for computing the derivative of C6 by finite difference
            double r = Coordination_Number_step[i];
            double t = Coordination_Number_step[j];

        //    inFile.open("UCHF_daug_C6.txt");
        //    inFile.open(UCHFPath.c_str());
            inFile.open("/home/chandemonium/MP2D_Program/Dispersion/UCHF_daug_C6.txt");
            string line;
            if (!inFile) {
            cerr << "Unable to open file UCHF_daug_C6.txt" << endl;;
            exit (1);
            }
            double GA_C6 =0;
            double sum_L_ij=0;
            double sum_numerator=0;
            // for C6 derivative
            double sum_C6ref_Lij =0;
            double sum_L_ij_deriv=0;
            double sum_C6ref_Lij_deriv=0;
            double sum_C6_deriv=0;

            // for finite difference 
            double sum_L_ij_step=0;
            double sum_C6ref_Lij_step=0;

            while (getline(inFile, line)) {
                istringstream ss(line);
                double var1, var2, var3, var4, var5;
                ss >> var1 >> var2 >> var3 >> var4 >> var5;
                if ((var1 == a && var2 == b && var3 == c && var4 == d) || (var1==b && var2==a && var3==d   && var4==c)) {
                //    cout << var5 << endl;
                UCHF[index] = var5;
                }
                if ((var1 == a && var2 == b) || (var1==b && var2==a)) {
              //  double L_ij = 
               //     cout << var8 << endl;
               //     cout << var9 << endl;
               //     cout << var10 << endl;
               //     double L_ij = exp(-4*(pow((e - round(var8)), 2) + pow((g - round(var9)), 2)));
                   
                    // unless the coordination numbers are equal, the contribution needs to be counted twice.
                    if (a==b && var3 != var4) {
                     //   L_ij = 2*exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9)))));    
                      //  L_ij = exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9))))) + exp(-4*(((e - round(var9))*(e - round(var9))) + ((g - round(var8))*(g - round(var8)))));
                        // The coordination numbers (hybridization numbers are not rounded to the nearest integer when calculating the Gaussian Averaged Dispersion Coefficients.
                        L_ij = exp(-4*(((e - var3)*(e - var3)) + ((g - var4)*(g - var4)))) + exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))));

                        // for finite difference
                        L_ij_step = exp(-4*(((r - var3)*(r - var3)) + ((t - var4)*(t - var4)))) + exp(-4*(((r - var4)*(r - var4)) + ((t - var3)*(t - var3))));


                 //*       L_ij_deriv = exp(-4*(((e - var3)*(e - var3)) + ((g - var4)*(g - var4))))*(-8*(((e-var3)*h)+((g-var4)*k))) + exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))))*(-8*(((e-var4)*k)+((g-var3)*h)));

                         L_ij_deriv = exp(-4*(((e - var3)*(e - var3)) + ((g - var4)*(g - var4))))*(-8*(((e-var3)*h)+((g-var4)*k))) + exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))))*(-8*(((e-var4)*h)+((g-var3)*k)));
                    } 
                    
                    // I have to use to else if statements to change the coordination numbers according to which atom comes first.
                    else if (var1 ==a && var2 ==b) {
                     //   L_ij = exp(-4*(((e - round(var8))*(e - round(var8))) + ((g - round(var9))*(g - round(var9)))));
                        L_ij = exp(-4*(((e - var3)*(e - var3)) + ((g - var4)*(g - var4))));

                        L_ij_step = exp(-4*(((r - var3)*(r - var3)) + ((t - var4)*(t - var4))));

                        L_ij_deriv = exp(-4*(((e - var3)*(e - var3)) + ((g - var4)*(g - var4))))*(-8*(((e-var3)*h)+((g-var4)*k)));
                    }
                    else if (var1==b && var2 ==a) {
                     //   L_ij = exp(-4*(((e - round(var9))*(e - round(var9))) + ((g - round(var8))*(g - round(var8)))));
                        L_ij = exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))));

                        L_ij_step = exp(-4*(((r - var4)*(r - var4)) + ((t - var3)*(t - var3))));

                  //*      L_ij_deriv = exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))))*(-8*(((e-var4)*k)+((g-var3)*h)));

                        L_ij_deriv = exp(-4*(((e - var4)*(e - var4)) + ((g - var3)*(g - var3))))*(-8*(((e-var4)*h)+((g-var3)*k)));
                    }
                    else {}

                    sum_L_ij += L_ij;
                    // for the C6 derivative
                    sum_C6ref_Lij += L_ij*var5;

                    sum_L_ij_deriv += L_ij_deriv;
                    sum_C6ref_Lij_deriv += L_ij_deriv*var5;

                    // for finite difference
                    sum_L_ij_step += L_ij_step;
                    sum_C6ref_Lij_step += L_ij_step*var5;

    //*                sum_C6_deriv += (sum_C6ref_Lij_deriv*sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv)/pow(sum_L_ij,2);

//                    cout << "sum L_ij = " << sum_L_ij << endl;
                //    cout << "L_ij = " << L_ij << endl;
                    double numerator = L_ij*var5;
                    sum_numerator += numerator;
//                    cout << "sum_numerator = " << sum_numerator << endl;
               //     cout << "e = " << e << endl;
                //    cout << "var8 = " << var8 << endl;
                 //   cout << "g = " << g << endl;
                 //   cout << "var9 = " << var9 << endl;
            // Neither equation for GA_C6 is correct as of 4/6/2018
                 //   GA_C6 += ((CKS[index] * L_ij) / L_ij);
                 //   GA_C6 += ((var10 * L_ij) / L_ij);
                 //   cout << "GA_C6 = " << GA_C6 << endl;
                GaussA_C6_UCHF[index] = sum_numerator/sum_L_ij;

      //*          C6_UCHF_deriv[index] = sum_C6_deriv;
                C6_UCHF_deriv[index] = (sum_C6ref_Lij_deriv*sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv)/pow(sum_L_ij,2);

                C6_UCHF_FD[index] = ((sum_C6ref_Lij_step/sum_L_ij_step) - (sum_numerator/sum_L_ij))/0.0000001;    

                }
            }
        //    cout << "C6 UCHF derivative " << C6_UCHF_deriv[index] << endl;
       //     cout << "The UCHF derivative by FD is " << C6_UCHF_FD[index] << endl;
            inFile.close(); 
        }       
    }    
    for (int i=0; i<combinations;i++) {
        
    //    cout << "UCHF = " << UCHF[i] << endl;
 //       cout << "UCHF = " << GaussA_C6_UCHF[i] << endl;
    }  
}

/*int factorial(int n) {
                    return (n==1 || n==0) ? 1 : factorial(n-1)*n;
} */


double Coord_Num::GetCKSEnergy() {

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 
    
    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();

    double *k;
    k = Dispersion::dispersion().CoordinateDerivs();

    
    double E6[combinations];
    d_E6_CKS = new double[combinations];
    double rcut = 0.72;
    double width = 0.20;
    double R_AB =0.0;
    double R_AB_deriv = 0.0;
    double R_AB_prime =0.0;
    double DampingTT = 0.0;
    double DampingTT_Deriv = 0.0;
    double E6_deriv_total = 0.0;
    E6_CKS_coord_derivs = new double[3*combinations];
    modR_AB = new double[combinations];
    double a_one = 0.9436334537945325;
    double a_two = 0.4802462930911932;
    double Total_E6_CKS = 0.0;
  //  double arg_1 = 0.0;
  //  double arg_2 = 0.0;
  //  double damp_sum = 0.0;
    int factorial = 0;
   // for (int i=0; i<combinations;i++) {
        int index = -1;
        int other_index = -1;
        for (int w=0; w<Ntot;w++) {
            for (int z=w+1; z<Ntot;z++) {
                index++;
                other_index +=3;
                double arg_1 = 0.0;
                double arg_2 = 0.0;
                double damp_sum = 0.0;
// Experimenting with implementing the derivatives for gradient calculations
                double arg_3 = 0.0;
                double arg_4 = 0.0;
                double damp_sum_deriv = 0.0;
         //       cout << *( *(p+w)+z) << endl;
         //       if ( (* ( *(p+w)+z)) <= RcovAB[w][z]*(rcut-(width/2))) {
                if ( (* ( *(p+w)+z)) <= r0AB[index]*(rcut-(width/2))) {
                 //   R_AB = rcut*RcovAB[w][z];
                //    cout << R_AB << endl;
                    R_AB = rcut*r0AB[index];
                    R_AB_deriv = 0.0;
                //    cout << R_AB << endl;
                }
          //      else if ( (* ( *(p+w)+z)) >= RcovAB[w][z]*(rcut+(width/2))) {
                else if ( (* ( *(p+w)+z)) >= r0AB[index]*(rcut+(width/2))) {
                    R_AB = (* ( *(p+w)+z));
                    R_AB_deriv = 1.0;
                //    cout << R_AB << endl;
                }
                else {
                    double r_prime = rcut*r0AB[index];
                    double w_prime = width*r0AB[index];
                    double function_x = ((* ( *(p+w)+z)) -(r_prime- (w_prime/2))) / (w_prime);
                //    cout << *(*(p+w)+z) << endl;
               //     cout << "function_x is" << function_x << endl;
                    double function_r = (-2.5*pow(function_x,8) + 10*pow(function_x,7) - 14*pow(function_x,6) + 7*pow(function_x,5))*w_prime;
                    double function_r_deriv = (-20*pow(function_x,7) + 70*pow(function_x,6) - 84*pow(function_x,5) + 35*pow(function_x,4)); //w_prime;
               //     cout << "function_r is " << function_r << endl;
               //     R_AB = rcut*RcovAB[w][z] + function_r;
                    R_AB = rcut*r0AB[index] + function_r;
                    R_AB_deriv = function_r_deriv;
                //   cout << R_AB << endl;
                } 
                // Store the modified values of R_AB
                modR_AB[index] = R_AB;
           
                // I stole this little function to do a factorial
          /*      int factorial(int n) {
                    return (n==1 || n==0) ? 1 : factorial(n-1)*n;
                } */
            
           //     cout << r0AB[index] << endl;                

                for (int b=0; b<=6; b++) {
                    if (b==0 || b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else {
                        factorial = 720;
                    }
                        
                  //  arg_1 = ((* ( *(p+w)+z))*a_one -a_two)*RcovAB[w][z];
           //         arg_1 = (R_AB*a_one - a_two)*RcovAB[w][z];
                 //   arg_1 = (R_AB*a_one - a_two)*r0AB[index];
                    arg_1 = ((r0AB[index]*a_one + a_two)*R_AB);
               //     cout << r0AB[index] << endl;
               //     cout << r0AB[index]*a_one -a_two << endl;
               //     cout << R_AB << endl;
              //      cout << r0AB[index] << endl;
              //      cout << "arg_1 " << arg_1 << endl;
                    arg_2 = pow(arg_1,b)/factorial;
             //       cout << "arg_2 " << arg_2 << endl;
                    damp_sum += arg_2;
         //           cout << "damp_sum " << damp_sum << endl;
                }
                DampingTT = 1 - exp(-1*arg_1)*damp_sum;
              //  cout << "arg_1 outside " << arg_1 << endl;
             //   cout << "the damping function is " << DampingTT << endl; 

// I need to calculate the derivative of the damping function for the gradients
                for (int b=1; b<=6; b++) {
                    if (b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else {
                        factorial = 720;
                    }
                        
  
                    arg_3 = ((r0AB[index]*a_one + a_two)*R_AB);
               
      //              arg_4 = pow(arg_3,b-1)/factorial;
                    arg_4 = (pow((r0AB[index]*a_one + a_two),b) * pow(R_AB,b-1))/factorial;
      //              cout << arg_4 << endl;
                    
             
                    damp_sum_deriv += arg_4*b;

// I need to calculate the derivative of C6 with respect to R
                
         
                }
// The original implementation of the damping function derivativie before I included the derivative of the double damping function.
  //              DampingTT_Deriv = (r0AB[index]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv);

                DampingTT_Deriv = ((r0AB[index]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv))*R_AB_deriv;
           //     cout << R_AB_deriv << endl;
        

 //       cout << damp_sum_deriv << endl;
//                cout << "DampingTT Derivative " << DampingTT_Deriv << endl;
                    
    
            

              // checking to make sure that I can print out the covalent radii from earlier  
              //    cout << RcovAB[w][z] << endl;
              //  cout << "C6_CKS = " << GaussA_C6_CKS[index] << endl;
              //  cout << "rAB = " << (* ( *(p+w)+z)) << endl;
                //double E_6 = (GaussA_C6_CKS[index] / pow((* ( *(p+w)+z)),6)) ;
                double E_6 =-1* ((GaussA_C6_CKS[index] / pow(R_AB*1.8897261254535,6)*627.5095));
           //    cout << GaussA_C6_CKS[index] << endl;
           //    cout << R_AB << endl;
           //     cout << R_AB*1.88973 << endl;
           //     cout << "the undamped energy is " << E_6 << endl;

// I might as well calculate the gradients while I'm at it
     //           double E6_deriv = (6*((GaussA_C6_CKS[index] / pow(R_AB*1.88973,7)*627.5095))*DampingTT) - ((GaussA_C6_CKS[index] / pow(R_AB*1.88973,6)*627.5095))*DampingTT_Deriv;
     //           cout << " E6_deriv " << E6_deriv << endl;
                double E6_deriv = ((6*GaussA_C6_CKS[index]*627.5095 / pow(R_AB,7))/pow(1.8897261254535,6))*R_AB_deriv;
         //       cout << E6_deriv << endl;
                double C6_Deriv = -1*((C6_CKS_deriv[index]/pow(R_AB*1.8897261254535,6)*627.5095));
                

          //      cout << "C6_Deriv " << C6_Deriv << endl;
// This was my expression before I included the derivative with respect to C6.
   //             double deriv6 = E_6*DampingTT_Deriv + E6_deriv*DampingTT;
         //       cout << " the derivative of the E6 CKS energy without C6 " << deriv6 << endl;

                double deriv6 = C6_Deriv*DampingTT + E_6*DampingTT_Deriv + E6_deriv*DampingTT;

                d_E6_CKS[index] = deriv6;
           //     d_E6_CKS[combinations + index ] = deriv6;
           //     d_E6_CKS[combinations+index] = deriv6;
           //     cout << "the derivaitve of the E6 CKS energy is " << deriv6 << endl;

           // Multiply dE/dR by dR/dx, dR/dy, and dR/dz to get dE/dxi, dE/dyi, dE/dzi.
                double deriv6_x = (deriv6* (*(k + (other_index-2))))*(1/R_AB);
           //     cout << "deriv6_x " << deriv6_x << endl;
                double deriv6_y = (deriv6* (*(k + (other_index-1))))*(1/R_AB);
            //    cout << "deriv6_y " << deriv6_y << endl;
                double deriv6_z = (deriv6* (*(k + (other_index))))*(1/R_AB);
            //    cout << "deriv6_z " << deriv6_z << endl;

                E6_CKS_coord_derivs[other_index-2] = deriv6_x;
                E6_CKS_coord_derivs[other_index-1] = deriv6_y;
                E6_CKS_coord_derivs[other_index] = deriv6_z;

         
        //        E6_deriv_total += deriv6;
        //        cout << "E6_deriv_total " << E6_deriv_total << endl;


                double E6_damped = E_6*DampingTT;
            //    cout << "the damped energy is " << E6_damped << endl;
            //    cout << E_6 << endl;
                E6[index] = E6_damped;
             //   cout << E6[index] << endl;
            //    E6[index] = E_6;
  //              Total_E6_CKS += E6_damped;
            //    cout << Total_E6 << endl;
                
            }
        }
        
        for (int i=0; i<combinations; i++) {
         //   cout << "E6 is " << E6[i] << endl;
            Total_E6_CKS += E6[i];
        }
        for (int i=0; i<combinations; i++) {
     //       cout << "E6 CKS coord derivs " << E6_CKS_coord_derivs[3*i] << " " << E6_CKS_coord_derivs[3*i+1] << " " << E6_CKS_coord_derivs[3*i+2] << endl;
       //     cout << "here goes nothin " << d_E6_CKS[i] << endl;
        }

  //      cout << "Total E6 CKS derivative " << E6_deriv_total << endl;
    //    cout << Total_E6_CKS << endl;
        return Total_E6_CKS;

}

double Coord_Num::GetCKSC8Energy() {

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    double *p;
    p = Dispersion::dispersion().GetMultipoleExpectationValue();
    // Test to make sure I'm printing out the r^2/r^4 expectation values.
   /* for (int i=0; i<Ntot; i++) {
        cout << *(p+i) << endl;
    }  */

    double *k;
    k = Dispersion::dispersion().CoordinateDerivs();

    int *h;
    h = Dispersion::dispersion().GetAtomicNumber();                              

 //   double C8CKS[combinations];                                    

  
   
//    double C8CKS[combinations];    
    C8CKS = new double[combinations];
    C8CKS_deriv = new double[combinations];
    int index = -1;
    int other_index = -1;
    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++) {
            index++;
            other_index +=3;
     //       cout << "CKS_C6 = " << GaussA_C6_CKS[index] << endl;
      //      cout << "the expectation value of atom A is " << *(p+i) << endl;
     //       cout << "the expectation value of atom B is " << *(p+j) << endl;
       //     cout << "the atomic number of atom A is " << *(h+i) << endl;
      //      cout << "the atomic number of atom B is " << *(h+j) << endl;
            // Grimme uses a constant s42 to scale Q_A and Q_B. I didn't see the value anywhere, so I used 
            // s42 = 0.5 which basically reproduces the Ruby4 values.
            double Q_A = 0.5*sqrt(*(h+i)) * *(p+i);
        //    cout << Q_A << endl;
            double Q_B = 0.5*sqrt(*(h+j)) * *(p+j);
        //    cout << Q_B << endl;
            double C8_CKS = 3*GaussA_C6_CKS[index]*sqrt(Q_A*Q_B);
       //     cout << C8_CKS << endl;

            C8CKS[index] = C8_CKS;
       //     cout << C8CKS[index] << endl;

             // Get the derivative of C8
            double C8_CKS_deriv = 3*C6_CKS_deriv[index]*sqrt(Q_A*Q_B);
            C8CKS_deriv[index] = C8_CKS_deriv;
    //        cout << "C8CKS derivative " << C8CKS_deriv[index] << endl;

            
        }
    }

    // I can get the C8s. Now I need to get the energies.
    double **l;
    l = Dispersion::dispersion().CalculateDistancesDifferently();

    double E8[combinations];
    d_E8_CKS = new double[combinations];
    double rcut = 0.72;
    double width = 0.20;
    double R_AB =0.0;
    double R_AB_deriv = 0.0;
    double DampingTT = 0.0;
    double DampingTT_Deriv = 0.0;
    double E8_deriv_total = 0.0;
    E8_CKS_coord_derivs = new double[3*combinations];
    double a_one = 0.9436334537945325;
    double a_two = 0.4802462930911932;
    double s_8 = 1.1873480299798238;
    double Total_E8_CKS = 0.0;
  //  double arg_1 = 0.0;
  //  double arg_2 = 0.0;
  //  double damp_sum = 0.0;
    int factorial = 0;
   // for (int i=0; i<combinations;i++) {
        int indexx = -1;
        int other_indexx = -1;
        for (int w=0; w<Ntot;w++) {
            for (int z=w+1; z<Ntot;z++) {
                indexx++;
                other_indexx +=3;
                double arg_1 = 0.0;
                double arg_2 = 0.0;
                double damp_sum = 0.0;
// I need to initiate some variables for the derivative
                double arg_3 = 0.0;
                double arg_4 = 0.0;
                double damp_sum_deriv = 0.0;
           //     if ( (* ( *(l+w)+z)) <= RcovAB[w][z]*(rcut-(width/2))) {
                if ( (* ( *(l+w)+z)) <= r0AB[indexx]*(rcut-(width/2))) {
                 //   R_AB = rcut*RcovAB[w][z];
                    R_AB = rcut*r0AB[indexx];
                    R_AB_deriv = 0.0;
                 //   cout << R_AB << endl;
                }
           //     else if ( (* ( *(l+w)+z)) >= RcovAB[w][z]*(rcut+(width/2))) {
                else if ( (* ( *(l+w)+z)) >= r0AB[indexx]*(rcut+(width/2))) {
                    R_AB = (* ( *(l+w)+z));
                    R_AB_deriv = 1.0;
                 //   cout << R_AB << endl;
                }
                else {
                    double r_prime = rcut*r0AB[indexx];
                    double w_prime = width*r0AB[indexx];
                    double function_x = ((* ( *(l+w)+z)) -(r_prime- (w_prime/2))) / (w_prime);
                    double function_r = (-2.5*pow(function_x,8) + 10*pow(function_x,7) - 14*pow(function_x,6) + 7*pow(function_x,5))*w_prime;
                    double function_r_deriv = (-20*pow(function_x,7) + 70*pow(function_x,6) - 84*pow(function_x,5) + 35*pow(function_x,4)); //w_prime;
                  //  R_AB = rcut*RcovAB[w][z] + function_r;
                    R_AB = rcut*r0AB[indexx] + function_r;
                    R_AB_deriv = function_r_deriv;
                 //   cout << R_AB << endl;
                } 
                // I stole this little function to do a factorial
          /*      int factorial(int n) {
                    return (n==1 || n==0) ? 1 : factorial(n-1)*n;
                } */
            
           //     cout << r0AB[index] << endl;                

                for (int b=0; b<=8; b++) {
                    if (b==0 || b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else if (b==6) {
                        factorial = 720;
                    }
                    else if (b==7) {
                        factorial = 5040;
                    }
                    else {
                        factorial = 40320;
                    }
                        
                  //  arg_1 = ((* ( *(p+w)+z))*a_one -a_two)*RcovAB[w][z];
           //         arg_1 = (R_AB*a_one - a_two)*RcovAB[w][z];
               //     arg_1 = (R_AB*a_one - a_two)*r0AB[index];
                    arg_1 = ((r0AB[indexx]*a_one + a_two)*R_AB);
                    arg_2 = pow(arg_1,b)/factorial;
                    damp_sum += arg_2;

                }
                DampingTT = 1 - exp(-1*arg_1)*damp_sum;
          //      cout << "the damping function is " << DampingTT << endl; 

// Here I'm implementing the first derivative of the E8 energy with respect to R_AB
                for (int b=1; b<=8; b++) {
                    if (b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else if (b==6) {
                        factorial = 720;
                    }
                    else if (b==7) {
                        factorial = 5040;
                    }
                    else {
                        factorial = 40320;
                    }
                        
  
                    arg_3 = ((r0AB[indexx]*a_one + a_two)*R_AB);
               
      //              arg_4 = pow(arg_3,b-1)/factorial;
                    arg_4 = (pow((r0AB[indexx]*a_one + a_two),b) * pow(R_AB,b-1))/factorial;
      //              cout << arg_4 << endl;
                    
             
                    damp_sum_deriv += arg_4*b;
         
                }
  // Before including the double damping derivative contribution
     //           DampingTT_Deriv = (r0AB[indexx]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv);

// with double damping derivative included
                DampingTT_Deriv = ((r0AB[indexx]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv))*R_AB_deriv;
         //       cout << damp_sum_deriv << endl;
            //    cout << "DampingTT Derivative " << DampingTT_Deriv << endl;
                    
                    
    
                double E_8 =-1*s_8*((C8CKS[indexx] / pow(R_AB*1.8897261254535,8))*627.5095);
    //            cout << "the undamped energy is " << E_8 << endl;
         //       cout << C8CKS[indexx] << endl;

//                double E8_deriv = (6*((C8CKS[indexx] / pow(R_AB*1.88973,7)*627.5095))*DampingTT) - ((C8CKS[indexx] / pow(R_AB*1.88973,6)*627.5095))*DampingTT_Deriv;
//                cout << " E8_deriv " << E8_deriv << endl;
                double E8_deriv = ((s_8*8*C8CKS[indexx]*627.5095 / pow(R_AB,9))/pow(1.8897261254535,8))*R_AB_deriv;
        //        cout << E8_deriv << endl;

                double C8_deriv = -1*s_8*((C8CKS_deriv[indexx] / pow(R_AB*1.8897261254535,8))*627.5095);

           // derivative without C8 derivative
 //               double deriv8 = E_8*DampingTT_Deriv + E8_deriv*DampingTT;
          //      cout << "CKS E8 derivative " << deriv8 << endl;
          // derivative of the CKS E8 energy with the derivative of C8
                double deriv8 = C8_deriv*DampingTT + E_8*DampingTT_Deriv + E8_deriv*DampingTT;
         //      cout << "the derivative of the E8 CKS energy is " << deriv8 << endl;
                d_E8_CKS[indexx] = deriv8;

                double deriv8_x = (deriv8* (*(k + (other_indexx-2))))/R_AB;
            //    cout << "deriv8_x " << deriv8_x << endl;
                double deriv8_y = (deriv8* (*(k + (other_indexx-1))))/R_AB;
            //    cout << "deriv8_y " << deriv8_y << endl;
                double deriv8_z = (deriv8* (*(k + (other_indexx))))/R_AB;
            //    cout << "deriv8_z " << deriv8_z << endl;

                E8_CKS_coord_derivs[other_indexx-2] = deriv8_x;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index-2] << endl;
                E8_CKS_coord_derivs[other_indexx-1] = deriv8_y;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index-1] << endl;
                E8_CKS_coord_derivs[other_indexx] = deriv8_z;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index] << endl;


         //       E8_deriv_total += E8_deriv;
                E8_deriv_total += deriv8;
//                cout << "E8_deriv_total " << E8_deriv_total << endl;

                double E8_damped = E_8*DampingTT;
     //           cout << "the damped energy is " << E8_damped << endl;
            //    cout << E_8 << endl;
                Total_E8_CKS += E8_damped;
               // cout << Total_E8_CKS << endl;
                
            }
        }
        for (int i=0; i<combinations; i++) {
 //           cout << "E8 CKS Coordinate derivatives " << E8_CKS_coord_derivs[3*i] << " "<< E8_CKS_coord_derivs[3*i+1] << " " << E8_CKS_coord_derivs[3*i+2] << endl;
      //      cout << "hoop hoop " << d_E8_CKS[i] << endl;
        }
  //      cout << " the derivative of E8 CKS is " << E8_deriv_total << endl;
        return Total_E8_CKS;
    

} 

double Coord_Num::GetUCHFEnergy() {

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 
    
    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();

    double *k;
    k = Dispersion::dispersion().CoordinateDerivs();

    
    double E6[combinations];
    d_E6_UCHF = new double[combinations];
    double rcut = 0.72;
    double width = 0.20;
    double R_AB =0.0;
    double R_AB_deriv = 0.0;
    double DampingTT = 0.0;
    double DampingTT_Deriv = 0.0;
    double E6_deriv_total = 0.0;
    E6_UCHF_coord_derivs = new double[3*combinations];
    double a_one = 0.9436334537945325;
    double a_two = 0.4802462930911932;
    double Total_E6_UCHF = 0.0;
  //  double arg_1 = 0.0;
  //  double arg_2 = 0.0;
  //  double damp_sum = 0.0;
    int factorial = 0;
   // for (int i=0; i<combinations;i++) {
        int index = -1;
        int other_index = -1;
        for (int w=0; w<Ntot;w++) {
            for (int z=w+1; z<Ntot;z++) {
                index++;
                other_index +=3;
                double arg_1 = 0.0;
                double arg_2 = 0.0;
                double damp_sum = 0.0;
                double arg_3 = 0.0;
                double arg_4 = 0.0;
                double damp_sum_deriv = 0.0;
            //    if ( (* ( *(p+w)+z)) <= RcovAB[w][z]*(rcut-(width/2))) {
                if ( (* ( *(p+w)+z)) <= r0AB[index]*(rcut-(width/2))) {
               //     R_AB = rcut*RcovAB[w][z];
                    R_AB = rcut*r0AB[index];
                    R_AB_deriv = 0.0;
                }
            //    else if ( (* ( *(p+w)+z)) >= RcovAB[w][z]*(rcut+(width/2))) {
                else if ( (* ( *(p+w)+z)) >= r0AB[index]*(rcut+(width/2))) {
                    R_AB = (* ( *(p+w)+z));
                    R_AB_deriv = 1.0;
                }
                else {
                    double r_prime = rcut*r0AB[index];
                    double w_prime = width*r0AB[index];
                    double function_x = ((* ( *(p+w)+z)) -(r_prime- (w_prime/2))) / (w_prime);
                    double function_r = (-2.5*pow(function_x,8) + 10*pow(function_x,7) - 14*pow(function_x,6) + 7*pow(function_x,5))*w_prime;
                    double function_r_deriv = (-20*pow(function_x,7) + 70*pow(function_x,6) - 84*pow(function_x,5) + 35*pow(function_x,4)); //w_prime;
                //    R_AB = rcut*RcovAB[w][z] + function_r;
                    R_AB = rcut*r0AB[index] + function_r;
                    R_AB_deriv = function_r_deriv;
                } 
                // I stole this little function to do a factorial
          /*      int factorial(int n) {
                    return (n==1 || n==0) ? 1 : factorial(n-1)*n;
                } */
            
           //     cout << r0AB[index] << endl;                

                for (int b=0; b<=6; b++) {
                    if (b==0 || b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else {
                        factorial = 720;
                    }
                        
                  //  arg_1 = ((* ( *(p+w)+z))*a_one -a_two)*RcovAB[w][z];
           //         arg_1 = (R_AB*a_one - a_two)*RcovAB[w][z];
               //     arg_1 = (R_AB*a_one - a_two)*r0AB[index];
                    arg_1 = ((r0AB[index]*a_one + a_two)*R_AB);
               //     cout << "arg_1 is " << arg_1 << endl;
                //    cout << "b is " << b << endl;
                    arg_2 = pow(arg_1,b)/factorial;
              //      cout << arg_2 << endl;
                    damp_sum += arg_2;
                }
      //          cout << "damp_sum is " << damp_sum << endl;
                DampingTT = 1 - exp(-1*arg_1)*damp_sum;
              //      cout << "arg_1 is " << arg_1 << endl;
      //          cout << "DampingTT is " << DampingTT << endl;
         //       cout << "the damping function is " << DampingTT << endl; 
              
                for (int b=1; b<=6; b++) {
                    if (b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else {
                        factorial = 720;
                    }
                        
  
                    arg_3 = ((r0AB[index]*a_one + a_two)*R_AB);
               
      //              arg_4 = pow(arg_3,b-1)/factorial;
                    arg_4 = (pow((r0AB[index]*a_one + a_two),b) * pow(R_AB,b-1))/factorial;
      //              cout << arg_4 << endl;
                    
             
                    damp_sum_deriv += arg_4*b;
         
                }

// Damping derivative before including the double damping function
     //           DampingTT_Deriv = (r0AB[index]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv);   

// Damping function derivative with double damping included
                                   
                DampingTT_Deriv = ((r0AB[index]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv))*R_AB_deriv; 
            

              // checking to make sure that I can print out the covalent radii from earlier  
              //    cout << RcovAB[w][z] << endl;
              //  cout << "C6_CKS = " << GaussA_C6_CKS[index] << endl;
              //  cout << "rAB = " << (* ( *(p+w)+z)) << endl;
                //double E_6 = (GaussA_C6_CKS[index] / pow((* ( *(p+w)+z)),6)) ;
                double E_6 =-1* ((GaussA_C6_UCHF[index] / pow(R_AB*1.8897261254535,6))*627.5095);
      //          cout << "the undamped energy is " << E_6 << endl;
      //          cout << "C6_UCHF " << GaussA_C6_UCHF[index] << endl;

                double E6_deriv = ((6*GaussA_C6_UCHF[index]*627.5095 / pow(R_AB,7))/pow(1.8897261254535,6))*R_AB_deriv;
         //       cout << E6_deriv << endl;

 // C6 derivative analytical
               double C6_Deriv = -1*((C6_UCHF_deriv[index]/pow(R_AB*1.8897261254535,6)*627.5095));
// C6 derivative finite difference
      //          double C6_Deriv = -1*((C6_UCHF_FD[index]/pow(R_AB*1.8897261254535,6)*627.5095));


          // Derivative before I included the derivative of C6
  //              double deriv6 = E_6*DampingTT_Deriv + E6_deriv*DampingTT;
          //      cout << "UCHF E6 derivative " << deriv6 << endl;

// The Derivative now includes the contribution from the derivative of the C6 coefficient.
                double deriv6 = C6_Deriv*DampingTT + E_6*DampingTT_Deriv + E6_deriv*DampingTT;
            //    cout << "the derivaitve of the E6 UCHF energy is " << deriv6 << endl;

                d_E6_UCHF[index] = deriv6;

                double deriv6_UCHF_x = (deriv6* (*(k + (other_index-2))))/R_AB;
           //     cout << "deriv6_UCHF_x " << deriv6_UCHF_x << endl;
                double deriv6_UCHF_y = (deriv6* (*(k + (other_index-1))))/R_AB;
           //     cout << "deriv6_UCHF_y " << deriv6_UCHF_y << endl;
                double deriv6_UCHF_z = (deriv6* (*(k + (other_index))))/R_AB;
           //     cout << "deriv6_UCHF_z " << deriv6_UCHF_z << endl;

                E6_UCHF_coord_derivs[other_index-2] = deriv6_UCHF_x;
                E6_UCHF_coord_derivs[other_index-1] = deriv6_UCHF_y;
                E6_UCHF_coord_derivs[other_index] = deriv6_UCHF_z;
                

                E6_deriv_total += deriv6;
//                cout << "E6_deriv_total " << E6_deriv_total << endl;


                double E6_damped = E_6*DampingTT;
     //           cout << "the damped energy is " << E6_damped << endl;
            //    cout << E_6 << endl;
                Total_E6_UCHF += E6_damped;
          //      cout << Total_E6 << endl;
                
            }
        }
        for (int i=0; i<combinations; i++) {
       //     cout << "E6 UCHF coord derivs " << E6_UCHF_coord_derivs[3*i] << " " << E6_UCHF_coord_derivs[3*i+1] << " " << E6_UCHF_coord_derivs[3*i+2] << endl;
        //   cout << "here ye here ye "<< d_E6_UCHF[i] << endl;
        }
     //   cout << "Total E6 UCHF deriviative " << E6_deriv_total << endl;
     //   cout << Total_E6_UCHF << endl;
        return Total_E6_UCHF;

}

double Coord_Num::GetUCHFC8Energy() {

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 

    double *p;
    p = Dispersion::dispersion().GetMultipoleExpectationValue();
    // Test to make sure I'm printing out the r^2/r^4 expectation values.
   /* for (int i=0; i<Ntot; i++) {
        cout << *(p+i) << endl;
    }  */

    double *k;
    k = Dispersion::dispersion().CoordinateDerivs();

    int *h;
    h = Dispersion::dispersion().GetAtomicNumber();                              

 //   double C8CKS[combinations];                                    

  
   
//    double C8CKS[combinations];    
    C8UCHF = new double[combinations];
    C8UCHF_deriv = new double[combinations];
    int index = -1;
    int other_index = -1;
    for (int i=0;i<Ntot;i++) {
        for (int j=i+1;j<Ntot;j++) {
            index++;
            other_index +=3;
     //       cout << "CKS_C6 = " << GaussA_C6_CKS[index] << endl;
      //      cout << "the expectation value of atom A is " << *(p+i) << endl;
     //       cout << "the expectation value of atom B is " << *(p+j) << endl;
       //     cout << "the atomic number of atom A is " << *(h+i) << endl;
      //      cout << "the atomic number of atom B is " << *(h+j) << endl;
            // Grimme uses a constant s42 to scale Q_A and Q_B. I didn't see the value anywhere, so I used 
            // s42 = 0.5 which basically reproduces the Ruby4 values.
            double Q_A = 0.5*sqrt(*(h+i)) * *(p+i);
        //    cout << Q_A << endl;
            double Q_B = 0.5*sqrt(*(h+j)) * *(p+j);
        //    cout << Q_B << endl;
            double C8_UCHF = 3*GaussA_C6_UCHF[index]*sqrt(Q_A*Q_B);
       //     cout << C8_CKS << endl;
            C8UCHF[index] = C8_UCHF;
       //     cout << C8CKS[index] << endl;

       // Get the derivative of C8
            double C8_UCHF_deriv = 3*C6_UCHF_deriv[index]*sqrt(Q_A*Q_B);
            C8UCHF_deriv[index] = C8_UCHF_deriv;

            
        }
    }

    // I can get the C8s. Now I need to get the energies.
    double **l;
    l = Dispersion::dispersion().CalculateDistancesDifferently();

    double E8[combinations];
    d_E8_UCHF = new double[combinations];
    double rcut = 0.72;
    double width = 0.20;
    double R_AB =0.0;
    double R_AB_deriv = 0.0;
    double DampingTT = 0.0;
    double DampingTT_Deriv = 0.0;
    double E8_deriv_total = 0.0;
    E8_UCHF_coord_derivs = new double[3*combinations];
    double a_one = 0.9436334537945325;
    double a_two = 0.4802462930911932;
    double s_8 = 1.1873480299798238;
    double Total_E8_UCHF = 0.0;
  //  double arg_1 = 0.0;
  //  double arg_2 = 0.0;
  //  double damp_sum = 0.0;
    int factorial = 0;
   // for (int i=0; i<combinations;i++) {
        int indexx = -1;
        int other_indexx = -1;
        for (int w=0; w<Ntot;w++) {
            for (int z=w+1; z<Ntot;z++) {
                indexx++;
                other_indexx +=3;
                double arg_1 = 0.0;
                double arg_2 = 0.0;
                double damp_sum = 0.0;
                double arg_3 = 0.0;
                double arg_4 = 0.0;
                double damp_sum_deriv = 0.0;
           //     if ( (* ( *(l+w)+z)) <= RcovAB[w][z]*(rcut-(width/2))) {
                if ( (* ( *(l+w)+z)) <= r0AB[indexx]*(rcut-(width/2))) {
                 //   R_AB = rcut*RcovAB[w][z];
                //    cout << R_AB << endl;
                    R_AB = rcut*r0AB[indexx];
                    R_AB_deriv = 0.0;
                 //   cout << r0AB[indexx] << endl;
                 //   cout << R_AB << endl;
                }
          //      else if ( (* ( *(l+w)+z)) >= RcovAB[w][z]*(rcut+(width/2))) {
                else if ( (* ( *(l+w)+z)) >= r0AB[indexx]*(rcut+(width/2))) {
                    R_AB = (* ( *(l+w)+z));
                    R_AB_deriv = 1.0;
                 //   cout << R_AB << endl;
                //    cout << r0AB[index] << endl;
                }
                else {
                    double r_prime = rcut*r0AB[indexx];
                    double w_prime = width*r0AB[indexx];
                    double function_x = ((* ( *(l+w)+z)) -(r_prime- (w_prime/2))) / (w_prime);
                    double function_r = (-2.5*pow(function_x,8) + 10*pow(function_x,7) - 14*pow(function_x,6) + 7*pow(function_x,5))*w_prime;
                    double function_r_deriv = (-20*pow(function_x,7) + 70*pow(function_x,6) - 84*pow(function_x,5) + 35*pow(function_x,4)); //w_prime;
              //      R_AB = rcut*RcovAB[w][z] + function_r;
                    R_AB = rcut*r0AB[indexx] + function_r;
                    R_AB_deriv = function_r_deriv;
                 //   cout << r0AB[index] << endl;
                //    cout << R_AB << endl;
                } 
                // I stole this little function to do a factorial
          /*      int factorial(int n) {
                    return (n==1 || n==0) ? 1 : factorial(n-1)*n;
                } */
            
           //     cout << r0AB[index] << endl;                

                for (int b=0; b<=8; b++) {
                    if (b==0 || b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else if (b==6) {
                        factorial = 720;
                    }
                    else if (b==7) {
                        factorial = 5040;
                    }
                    else {
                        factorial = 40320;
                    }
                        
                  //  arg_1 = ((* ( *(p+w)+z))*a_one -a_two)*RcovAB[w][z];
           //         arg_1 = (R_AB*a_one - a_two)*RcovAB[w][z];
               //     arg_1 = (R_AB*a_one - a_two)*r0AB[index];
                    arg_1 = ((r0AB[indexx]*a_one + a_two)*R_AB);
                    arg_2 = pow(arg_1,b)/factorial;
                    damp_sum += arg_2;
                }
                DampingTT = 1 - exp(-1*arg_1)*damp_sum;

                for (int b=1; b<=8; b++) {
                    if (b==1) {
                        factorial = 1;
                    }
                    else if (b==2) {
                        factorial = 2;
                    }
                    else if (b==3) {
                        factorial = 6;
                    }
                    else if (b==4) {
                        factorial = 24;
                    }
                    else if (b==5) {
                        factorial = 120;
                    }
                    else if (b==6) {
                        factorial = 720;
                    }
                    else if (b==7) {
                        factorial = 5040;
                    }
                    else {
                        factorial = 40320;
                    }
                        
  
                    arg_3 = ((r0AB[indexx]*a_one + a_two)*R_AB);
               
      //              arg_4 = pow(arg_3,b-1)/factorial;
                    arg_4 = (pow((r0AB[indexx]*a_one + a_two),b) * pow(R_AB,b-1))/factorial;
      //              cout << arg_4 << endl;
          //          cout << R_AB << endl;
                    
             
                    damp_sum_deriv += arg_4*b;
         
                }
// Damping function derivative without the double damping function contribution
    //            DampingTT_Deriv = (r0AB[indexx]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv);

// Damping function derivative with the double damping function contribution
                DampingTT_Deriv = ((r0AB[indexx]*a_one + a_two)*exp(-1*arg_3)*damp_sum - (exp(-1*arg_3)*damp_sum_deriv))*R_AB_deriv;  
              //  cout << R_AB_deriv << endl;
         


        //       cout << damp_sum_deriv << endl;
            //    cout << "DampingTT Derivative " << DampingTT_Deriv << endl;
                    
                                    


                 //   cout << arg_1 << endl;
                 //   cout << damp_sum << endl;
          //      cout << "the damping function is " << DampingTT << endl; 
                    
    
                double E_8 =-1*s_8*((C8UCHF[indexx] / pow(R_AB*1.8897261254535,8))*627.5095);
    //              cout << C8UCHF[indexx] << endl;
    //            cout << "the undamped energy is " << E_8 << endl;

                double E8_deriv = ((s_8*8*C8UCHF[indexx]*627.5095 / pow(R_AB,9))/pow(1.8897261254535,8))*R_AB_deriv;
        //        cout << E8_deriv << endl;

                double C8_deriv = -1*s_8*((C8UCHF_deriv[indexx] / pow(R_AB*1.8897261254535,8))*627.5095);

        //      derivative before I included the C8 derivative
   //             double deriv8 = E_8*DampingTT_Deriv + E8_deriv*DampingTT;
            //    cout << " E8 UCHF derivative " << deriv8 << endl;
       // derivative of the CKS E8 energy with the derivative of C8
                double deriv8 = C8_deriv*DampingTT + E_8*DampingTT_Deriv + E8_deriv*DampingTT;
            //    cout << "the derivative of the E8 UCHF energy is " << deriv8 << endl;

                d_E8_UCHF[indexx] = deriv8;

                double deriv8_UCHF_x = (deriv8* (*(k + (other_indexx-2))))/R_AB;
             //   cout << "deriv8_UCHF_x " << deriv8_UCHF_x << endl;
                double deriv8_UCHF_y = (deriv8* (*(k + (other_indexx-1))))/R_AB;
            //    cout << "deriv8_UCHF_y " << deriv8_UCHF_y << endl;
                double deriv8_UCHF_z = (deriv8* (*(k + (other_indexx))))/R_AB;
             //   cout << "deriv8_UCHF_z " << deriv8_UCHF_z << endl;

                E8_UCHF_coord_derivs[other_indexx-2] = deriv8_UCHF_x;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index-2] << endl;
                E8_UCHF_coord_derivs[other_indexx-1] = deriv8_UCHF_y;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index-1] << endl;
                E8_UCHF_coord_derivs[other_indexx] = deriv8_UCHF_z;
             //   cout << "works" << E8_UCHF_coord_derivs[other_index] << endl;


                E8_deriv_total += deriv8;
//                cout << "E8_deriv_total " << E8_deriv_total << endl;

                double E8_damped = E_8*DampingTT;
     //           cout << "the damped energy is " << E8_damped << endl;
            //    cout << E_8 << endl;
                Total_E8_UCHF += E8_damped;
          //      cout << Total_E8_UCHF << endl;
                
            }
         
        }
        for (int i=0; i<combinations; i++) {
     //       cout << "E8 UCHF coord derivs " << E8_UCHF_coord_derivs[3*i] << " " << E8_UCHF_coord_derivs[3*i+1] << " " << E8_UCHF_coord_derivs[3*i+2] << endl;
       //     cout << "Lebron is the GOAT " << d_E8_UCHF[i] << endl;
        }
     
   //     cout << " the derivative of E8 UCHF is " << E8_deriv_total << endl;
        return Total_E8_UCHF;
    

} 

double Coord_Num::MP2DDispersionCorrection() {
    double E_CKS = GetCKSEnergy();
   //     cout << E_CKS << endl;
    double E_UCHF = GetUCHFEnergy();  
   //     cout << E_UCHF << endl;
    double E8_CKS = GetCKSC8Energy();
   //     cout << E8_CKS << endl;
    double E8_UCHF = GetUCHFC8Energy();
   //     cout << E8_UCHF << endl;
    double Result = (E_CKS + E8_CKS) - (E_UCHF + E8_UCHF);
    cout << "Edisp kcal/mol " << Result << endl;  
  //  cout << Result << endl;
    cout << "normal termination of mp2d" << endl;

//Also I am going to suppress this for the now for the sake of the Psi4 integration 5/29/2018
 /*   ofstream myfile;
    myfile.open ("Dispersion_Energy.txt");
    myfile << "Dispersion Energy (kcal/mol) = " << Result << endl;
    myfile << "MP2D Terminated Successfully \n" << endl; 
    myfile.close(); */
   
}

double Coord_Num::GetGradient() {

    int Ntot = Dispersion::dispersion().GetTotalNumberOfAtoms();

    int combinations = 0;
    for (int i=0;i<Ntot; i++) {
        combinations += i;
    } 
    
    double **p;
    p = Dispersion::dispersion().CalculateDistancesDifferently();

    double **k;
    k = Dispersion::dispersion().CoordinateDerivs_2();

    string *b;
    b = Dispersion::dispersion().GetSymbols();

    double rcut = 0.72;
    double width = 0.20;
    double R_AB =0.0;
    int index = -1;
    int other_index = -1;

 //   double Distances[combinations];
/*    double Distances[Ntot*Ntot];
    double New_r0AB[2*combinations];

    
    for (int i=0; i<combinations; i++) {
        New_r0AB[i] = r0AB[i];
        New_r0AB[combinations + i] = r0AB[i];
    }
    
    for (int w=0; w<Ntot;w++) {
          //  for (int z=w+1; z<Ntot;z++) {
        for (int z=0; z<Ntot; z++) {
           //  index++;
            // cout << index << endl;
            //    index++;
            if (w==z) {
             //   break;
                R_AB = 0;
            }
            else if (w!=z) {  
                index++; 
             //   cout << *(*(p+w)+z) << endl;
         //      cout << index << endl;
            //    other_index +=3;
            //    if (w < z || w > z) {
            //        index++;
            //    }
            //    else {}
            //    if ( (* ( *(p+w)+z)) <= RcovAB[w][z]*(rcut-(width/2))) {
             //   if (w==z) {
             //       R_AB = 0;
             //   }
             //   else if ( (* ( *(p+w)+z)) <= r0AB[index]*(rcut-(width/2))) {
             //   else if ( (* ( *(p+w)+z)) <= New_r0AB[index]*(rcut-(width/2))) {
                if ( (* ( *(p+w)+z)) <= New_r0AB[index]*(rcut-(width/2))) {
                  //  cout << (* ( *(p+w)+z)) << " " << New_r0AB[index]*(rcut-(width/2)) << endl;
                  //  cout << New_r0AB[index] << endl;
             //       cout << *(*(p+w)+z) << endl;
                      
               //     R_AB = rcut*RcovAB[w][z];
           //         R_AB = rcut*r0AB[index];
                    R_AB = rcut*New_r0AB[index];
                    Distances[index] = R_AB;
                //    cout << "R_AB" << R_AB << endl;
                }
            //    else if ( (* ( *(p+w)+z)) >= RcovAB[w][z]*(rcut+(width/2))) {
          //      else if ( (* ( *(p+w)+z)) >= r0AB[index]*(rcut+(width/2))) {
                else if ( (* ( *(p+w)+z)) >= New_r0AB[index]*(rcut+(width/2))) {
                 //   cout << New_r0AB[index] << endl;
               //     cout << *(*(p+w)+z) << endl;
                    
                    R_AB = (* ( *(p+w)+z));
                    Distances[index] = R_AB;
               //     cout << "R_AB " << R_AB << endl;
                }
                else {
                //    cout << New_r0AB[index] << endl;
              //      cout << *(*(p+w)+z) << endl;
                   
                    
             //       double r_prime = rcut*r0AB[index];
                    double r_prime = rcut*New_r0AB[index];
             //       double w_prime = width*r0AB[index];
                    double w_prime = width*New_r0AB[index];
                    double function_x = ((* ( *(p+w)+z)) -(r_prime- (w_prime/2))) / (w_prime);
                    double function_r = (-2.5*pow(function_x,8) + 10*pow(function_x,7) - 14*pow(function_x,6) + 7*pow(function_x,5))*w_prime;
                //    R_AB = rcut*RcovAB[w][z] + function_r;
                    R_AB = rcut*New_r0AB[index] + function_r;
                    Distances[index] = R_AB;
        //            R_AB = rcut*r0AB[index] + function_r;
                 //   cout << " R_AB " << R_AB << endl;
                } 
           //     Distances[index] = R_AB;
          //      cout << "Distances " <<Distances[index] << endl;
        //        cout << "cutoff radii " << New_r0AB[index] << endl;
               // cout << "cutoff distances " << New_r0AB[index] << endl;
               // cout << "RAB " << *(*(p+w)+z) << endl;
            }
            else {}
            cout << "RAB " << *(*(p+w)+z) << endl;
        //    cout << index << endl;
        }
    }*/
 /*   for (int i =0; i<combinations; i++ ) {

        cout << "E6 deriv " << d_E6_CKS[i] << endl;
    }*/


    double new_d_E6_CKS[Ntot][Ntot];
    double new_d_E8_CKS[Ntot][Ntot];
    double new_d_E6_UCHF[Ntot][Ntot];
    double new_d_E8_UCHF[Ntot][Ntot];
    
    double Distances[Ntot][Ntot];
    int l =-1;

    for (int w=0; w<Ntot;w++) {
          //  for (int z=w+1; z<Ntot;z++) {
        for (int z=0; z<Ntot; z++) {
            if (w == z ) {
                new_d_E6_CKS[w][z] = 0.0;
                new_d_E8_CKS[w][z] = 0.0;
                new_d_E6_UCHF[w][z] = 0.0;
                new_d_E8_UCHF[w][z] = 0.0;
                Distances[w][z] = 0.0;
            }
            else if (w < z) {
                l++;
                new_d_E6_CKS[w][z] = d_E6_CKS[l];
                new_d_E8_CKS[w][z] = d_E8_CKS[l];
                new_d_E6_UCHF[w][z] = d_E6_UCHF[l];
                new_d_E8_UCHF[w][z] = d_E8_UCHF[l];
                Distances[w][z] = modR_AB[l];
            }
            else if (w > z) {
                new_d_E6_CKS[w][z] = new_d_E6_CKS[z][w];
                new_d_E8_CKS[w][z] = new_d_E8_CKS[z][w];
                new_d_E6_UCHF[w][z] = new_d_E6_UCHF[z][w];
                new_d_E8_UCHF[w][z] = new_d_E8_UCHF[z][w];
                Distances[w][z] = Distances[z][w];
            }
            else {}
         //       cout << "d_E6_CKS " << new_d_E6_CKS[w][z] << endl;
        //    cout << "Distances " << Distances[w][z] << endl;

        }
    }

    vector<double> Gradient_V(3*Ntot);
    
    double Gradient_X[Ntot];
    double Gradient_Y[Ntot];
    double Gradient_Z[Ntot];



  //  int new_index = -1;
    double deriv6_CKS_x =0;
    double deriv6_CKS_y =0;
    double deriv6_CKS_z =0;

    double deriv8_CKS_x =0;
    double deriv8_CKS_y =0;
    double deriv8_CKS_z =0;

    double deriv6_UCHF_x =0;
    double deriv6_UCHF_y =0;
    double deriv6_UCHF_z =0;

    double deriv8_UCHF_x =0;
    double deriv8_UCHF_y =0;
    double deriv8_UCHF_z =0;

    double CKS_Gradient[Ntot][3*Ntot];
    double UCHF_Gradient[Ntot][3*Ntot];



    for (int w=0; w<Ntot; w++) {
        Gradient_X[w] = 0.0;
        Gradient_Y[w] = 0.0;
        Gradient_Z[w] = 0.0;
        for (int z=0; z<Ntot; z++) {
           //     new_index++;
                if (w == z) {
                    
                    deriv6_CKS_x =0;
                    deriv6_CKS_y =0;
                    deriv6_CKS_z =0;

                    deriv8_CKS_x =0;
                    deriv8_CKS_y =0;
                    deriv8_CKS_z =0;

                    deriv8_UCHF_x =0;
                    deriv8_UCHF_x =0;
                    deriv8_UCHF_x =0;

                    CKS_Gradient[w][3*z] = 0.0;
                    CKS_Gradient[w][3*z+1] = 0.0;
                    CKS_Gradient[w][3*z+2] = 0.0;
                    UCHF_Gradient[w][3*z] = 0.0;
                    UCHF_Gradient[w][3*z+1] = 0.0;
                    UCHF_Gradient[w][3*z+2] = 0.0;

                    
                }
                else {
                //    new_index++;
              //      cout << new_index << endl;
                    
                
               //     double deriv6_CKS_x = (d_E6_CKS[index]* (*(*(k + (other_index-2))+z)))/R_AB;
                    deriv6_CKS_x = (new_d_E6_CKS[w][z]* (*(*(k + w)+(3*z))))/Distances[w][z];
                    deriv8_CKS_x = (new_d_E8_CKS[w][z]* (*(*(k + w)+(3*z))))/Distances[w][z];
                    deriv6_UCHF_x = (new_d_E6_UCHF[w][z]* (*(*(k + w)+(3*z))))/Distances[w][z];
                    deriv8_UCHF_x = (new_d_E8_UCHF[w][z]* (*(*(k + w)+(3*z))))/Distances[w][z];
                    CKS_Gradient[w][3*z] = deriv6_CKS_x + deriv8_CKS_x;
                    UCHF_Gradient[w][3*z] = deriv6_UCHF_x + deriv8_UCHF_x;
          //          cout << "d_E6_CKS " << new_d_E6_CKS[w][z] << endl;
          //          cout << "d_x " << *(*(k + w) +(3*z)) << endl;
          //          cout << "the distance " << Distances[w][z] << endl;
               //     cout << *(*(k + w)+(3*z)) << endl;
                //    cout << (*(*(k + w)+z+2)) << endl;
                //    cout << "Distance " << Distances[w][z] << endl;
       //             cout << "deriv6_CKS_x " << deriv6_CKS_x << endl;
                //    cout << Distances[new_index] << endl;
                //    double deriv6_CKS_y = (d_E6_CKS[index]* (*(*(k + (other_index-1))+z)))/R_AB;
                    deriv6_CKS_y = (new_d_E6_CKS[w][z]* (*(*(k + w)+(3*z+1))))/Distances[w][z];
                    deriv8_CKS_y = (new_d_E8_CKS[w][z]* (*(*(k + w)+(3*z+1))))/Distances[w][z];
                    deriv6_UCHF_y = (new_d_E6_UCHF[w][z]* (*(*(k + w)+(3*z+1))))/Distances[w][z];
                    deriv8_UCHF_y = (new_d_E8_UCHF[w][z]* (*(*(k + w)+(3*z+1))))/Distances[w][z];
                    CKS_Gradient[w][3*z+1] = deriv6_CKS_y + deriv8_CKS_y;
                    UCHF_Gradient[w][3*z+1] = deriv6_UCHF_y + deriv8_UCHF_y;

            //        cout << "d_E8_CKS " << new_d_E8_CKS[w][z] << endl;
            //        cout << "d_x " << *(*(k + w) +(3*z)) << endl;
            //        cout << "the distance " << Distances[w][z] << endl;
                //    cout << *(*(k + w)+(3*z+1)) << endl;
       //             cout << "deriv6_CKS_y " << deriv6_CKS_y << endl;
               //     cout << Distances[new_index] << endl;
              //      double deriv6_CKS_z = (d_E6_CKS[index]* (*(*(k + (other_index))+z)))/R_AB;
                    deriv6_CKS_z = (new_d_E6_CKS[w][z]* (*(*(k + w)+(3*z+2))))/Distances[w][z];
                    deriv8_CKS_z = (new_d_E8_CKS[w][z]* (*(*(k + w)+(3*z+2))))/Distances[w][z];
                    deriv6_UCHF_z = (new_d_E6_UCHF[w][z]* (*(*(k + w)+(3*z+2))))/Distances[w][z];
                    deriv8_UCHF_z = (new_d_E8_UCHF[w][z]* (*(*(k + w)+(3*z+2))))/Distances[w][z];
                    CKS_Gradient[w][3*z+2] = deriv6_CKS_z + deriv8_CKS_z;
                    UCHF_Gradient[w][3*z+2] = deriv6_UCHF_z + deriv8_UCHF_z;
                //    cout << *(*(k + w)+(3*z+2)) << endl;
       //             cout << "deriv6_CKS_z " << deriv6_CKS_z << endl;
              //      cout << Distances[new_index] << endl;

                    double argument_x = (deriv6_CKS_x + deriv8_CKS_x) - (deriv6_UCHF_x + deriv8_UCHF_x); 
                 //   cout << "arugment x " << argument_x << endl;
                    double argument_y = (deriv6_CKS_y + deriv8_CKS_y) - (deriv6_UCHF_y + deriv8_UCHF_y);
                //    cout << "arugment y " << argument_y << endl;
                    double argument_z = (deriv6_CKS_z + deriv8_CKS_z) - (deriv6_UCHF_z + deriv8_UCHF_z);
                //    cout << "arugment z " << argument_z << endl;
                    Gradient_X[w] += argument_x;
                    Gradient_Y[w] += argument_y;
                    Gradient_Z[w] += argument_z;

                    Gradient_V[3*w] += argument_x;
                    Gradient_V[3*w+1] += argument_y;
                    Gradient_V[3*w+2] += argument_z;
                    
                    /// Everything works up until this point, but I'm going to try and save my gradient as a vector
                    
                   
                }
                
               
             //   cout << " Coordinate CKS " << CKS_Gradient[w][3*z] << " " << CKS_Gradient[w][3*z+1] << " " << CKS_Gradient[w][3*z+2] << endl;
             //   cout << " Coordinate UCHF " << UCHF_Gradient[w][3*z] << " " << UCHF_Gradient[w][3*z+1] << " " << UCHF_Gradient[w][3*z+2] << endl;
                  
        }
    } 

    // For now I am going to suppress this for the sake of the psi4 integration process 5/29/2018

    ofstream myfile;
    myfile.open ("mp2d_gradient");
//    myfile << "Dispersion Energy (kcal/mol) = " << Result << endl;
//    myfile << "MP2D Terminated Successfully \n" << endl; 
//    myfile << "Gradient (kcal/mol/Angstrom): \n" << endl;
//    myfile.close();
    double sum_sq = 0.0;
//    printf ("Gradient (kcal/mol/Angstrom): \n");
    for (int i =0; i<Ntot; i++ ) {
      //  cout << Gradient_X[i] << "    " << Gradient_Y[i] << "    " << Gradient_Z[i] << endl;
     //   cout << Gradient_V[3*i] << "    " << Gradient_V[3*i+1] << "    " << Gradient_V[3*i+2] << endl;
        
        printf ("%*f %*f %*f \n", 15, Gradient_V[3*i]/(627.5095*1.8897259886), 15, Gradient_V[3*i+1]/(627.5095*1.8897259886), 15, Gradient_V[3*i+2]/(627.5095*1.8897259886));
    //    cout << *(b + i) << endl;
//        myfile << *(b +i) << "     " << Gradient_V[3*i] << "     " << Gradient_V[3*i+1] << "     " << Gradient_V[3*i+2] << endl;
        // The gradients are converted to hartree/Bohr so that PSI4 can use the gradients
        myfile << Gradient_V[3*i]/(627.5095*1.8897259886) << "     " << Gradient_V[3*i+1]/(627.5095*1.8897259886) << "     " << Gradient_V[3*i+2]/(627.5095*1.8897259886) << endl;
       
        
        sum_sq += (Gradient_V[3*i]*Gradient_V[3*i] + Gradient_V[3*i+1]*Gradient_V[3*i+1] + Gradient_V[3*i+2]*Gradient_V[3*i+2]);
        
    }
    double RMS = sqrt(sum_sq/(3*Ntot));
  //  printf (" \nRMS Gradient (kcal/mol/Angstrom): %f \n\n", RMS); 
  //  myfile << "\nRMS Gradient (kcal/mol/Angstrom): \n\n" << RMS << endl;
    myfile.close();
    

}












