#ifndef COMMONPHYSICS_H
#define COMMONPHYSICS_H

#include "commonConstants.h"

void mass_to_mole_frac(double* X,           // Mass fractions to be altered
                       const double* Y,     // Mole fractions input
                       const double* MW,    // Molecular weights 
                       const int& NS) {      // Number of species

    double denom = 0.0;
    for (int i = 0; i < NS; ++i) 
        denom += Y[i] / MW[i];    

    for (int i = 0; i < NS; ++i) 
        X[i] = (Y[i] / MW[i]) / denom;
    
}

void mixture_gas_constant(double& R_mix,        // Mixture gas constant to be altered
                          const double* Y,      // Mass fractions input
                          const double* MW,     // Molecular weights
                          const int& NS) {      // Number of species

    R_mix = 0.0;
    for (int i = 0; i < NS; ++i) 
        R_mix += Y[i] * (ugconn / MW[i]);
    
}

void mole_to_mass_frac(double* Y,           // Mass fractions to be altered
                       const double* X,     // Mole fractions input
                       const double* MW,    // Molecular weights
                       const int& NS) {     // Number of species

    double denom = 0.0;
    for (int i = 0; i < NS; ++i) 
        denom += X[i] * MW[i];
    

    for (int i = 0; i < NS; ++i) 
        Y[i] = (X[i] * MW[i]) / denom;    
}

void total_density(double& rho_tot,         // Total density to be altered
                   const double* rho_s,
                   const int& NS) {     // Mass fractions input

    rho_tot = 0.0;
    for (int i = 0; i < NS; ++i) 
        rho_tot += rho_s[i];

}

void compute_pressure(double& p,            // Pressure to be altered
                      const double& rho,    // Density
                      const double& R_mix,  // Mixture gas constant
                      const double& T) {    // Temperature 

    p = rho * R_mix * T;
}

#endif