#ifndef ENERGY_H
#define ENERGY_H

#include "dataStructures.h"
#include "common_functions.h"
#include "common.h"

void initialize_energy(ReactionSet& RS, double& e_total, double& e_vib, double* rho_s, double* Ts) {

    e_total = 0.0;
    e_vib = 0.0;

    for (int i = 0; i < RS.n_species; ++i) {

        // Add translational energy
        e_total += rho_s[i] * 3.0 / 2.0 * RS.species[i].R * Ts[0];

        // Add rotational energy if molecule
        if (RS.species[i].mol) {
            e_total += rho_s[i] * RS.species[i].R * Ts[0];
        }

        // Add vibrational energy if molecule
        if (RS.species[i].mol) {
            double theta_v = RS.species[i].theta_v;
            double ev = RS.species[i].R * theta_v / (exp(theta_v / Ts[1]) - 1.0);
            e_total += rho_s[i] * ev;
            e_vib += rho_s[i] * ev;
        }

        // Add formation energy 
        e_total += rho_s[i] * RS.species[i].ef;
    }
}




#endif