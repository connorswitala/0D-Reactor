#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "common.h"


// Linear interpolation function
double lerp(double y0, double x0, double y1, double x1, double x) {
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

void mass_to_mole(const ReactionSet& RS, const double* Y, double* X) {

    double denom = 0.0;
    for (int i = 0; i < RS.n_species; ++i) {
        denom += Y[i] / RS.species[i].mw;
    }

    for (int i = 0; i < RS.n_species; ++i) {
        X[i] = (Y[i] / RS.species[i].mw) / denom;
    }
}

void mixture_R(const ReactionSet& RS, const double* Y, double& R_mix) {

    R_mix = 0.0;
    for (int i = 0; i < RS.n_species; ++i) {
        R_mix += Y[i] * (ugconn / RS.species[i].mw);
    }
}

void mole_to_mass(const ReactionSet& RS, const double* X, double* Y) {

    double denom = 0.0;
    for (int i = 0; i < RS.n_species; ++i) {
        denom += X[i] * RS.species[i].mw;
    }

    for (int i = 0; i < RS.n_species; ++i) {
        Y[i] = (X[i] * RS.species[i].mw) / denom;
    }
}

void total_density(const ReactionSet& RS, const double* Y, double& rho_tot) {

    rho_tot = 0.0;
    for (int i = 0; i < RS.n_species; ++i) {
        rho_tot += Y[i];
    }
}

void compute_pressure(double& p, const double& rho, const double& R_mix, const double& T) {
    p = rho * R_mix * T;
}


#endif