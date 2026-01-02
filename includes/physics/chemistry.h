#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "../xml/reader.h"
#include "dataStructures.h"

/**
 *      Functions:
 *      
 *      find_nd_level() - find number density level for Park Keq coefficients
 *      
 *      compute_rates() - compute chemical production rates
 */


int find_nd_level(const std::vector<double>& Ns, 
                  double n) {

    const int N = (int)Ns.size();
    if (N < 2) return 0;

    if (n <= Ns[0]) return 0;
    if (n >= Ns[N-1]) return N-2;

    // find j with Ns[j] <= n < Ns[j+1]
    int j = 0;
    while (j+1 < N && !(Ns[j] <= n && n < Ns[j+1])) ++j;
    if (j > N-2) j = N-2;
    return j;
}


void compute_rates(double* rates,   // Production rates array
    ReactionSet& RS,                // Reaction set
    double* rho_s,                  // Species densities
    double* Ts,                     // Temperatures
    double& P) {                    // Pressure in Pascals
        

    // Molar concentration array
    std::vector<double> C(RS.n_species);

    // Compute concentrations
    for (int i = 0; i < RS.n_species; ++i) 
        C[i] = (rho_s[i] / RS.MWs[i]) * 1e-6; // mol/cm^3

    // Calculate rates from each reaction
    for (auto& r : RS.reactions) {

        double T, gamma, n;     // Intermediate variables for computation

        std::vector<double> A(r.Keq_N);         // Vector for linear interpolated Keq
        n = P / (boltzmann * Ts[0]) * 1e-6;     // Compute number density for curve fits
        int level = find_nd_level(r.Ns, n);     // Find number density level
 
        // Find Keq coefficients by linear interpolation
        for (int i = 0; i < r.Keq_N; ++i) {
            A[i] = lerp(r.Keq_As[level * r.Keq_N + i], 
                r.Ns[level], 
                r.Keq_As[(level + 1) * r.Keq_N + i], 
                r.Ns[level + 1], 
                n);
        }

        T = pow(Ts[0], r.Texp) * pow(Ts[1], (1.0 - r.Texp)); // Reaction temperature

        double exp_term = A[0] * (T / 10000.0)
                            + A[1]
                            + A[2] * log(10000.0 / T)
                            + A[3] * (10000.0 / T)
                            + A[4] * (10000.0 * 10000.0) / (T * T);

        double Keq = exp(exp_term); // Equilibrium constant

        double R;           // Reaction rate without (nu_b - nu_f) in front
        double f = 1.0;     // Forward part of reaction
        double b = 1.0;     // Backward part of reaction
        double nu;          // Stoichiometric coefficient

        gamma = r.C * pow(T, r.N) * exp(-r.Ea/T); // Base forward reaction coefficient - k_f = cm^3 / [mol -s]
            
        // Compute forward part
        for (int i = 0; i < RS.n_species; ++i) {
            nu = RS.nus_f[RS.n_species * r.id + i];
            if (nu == 0)
                continue;
            f *= pow(C[i], nu);
        }

        // Compute backward part
        for (int i = 0; i < RS.n_species; ++i) {
            nu = RS.nus_b[RS.n_species * r.id + i];
            if (nu == 0)
                continue;
            b *= pow(C[i], nu);
        }

        // Compute reaction rate depending on reaction type
        if (r.third_body) {
            double M_eff = 0.0;
            for (int i = 0; i < RS.n_species; ++i) {
                M_eff += C[i] * r.efficiencies[i];
            }

            R = gamma * M_eff * (f - b / Keq);
        }
        else {
            R = gamma * (f - b / Keq); // mol / (cm^3 / s)
        }    

        // Compute dXdt 
        for (int i = 0; i < RS.n_species; ++i) {
            int nu_f = RS.nus_f[RS.n_species * r.id + i];
            int nu_b = RS.nus_b[RS.n_species * r.id + i];

            r.dXdt[i] = R * (nu_b - nu_f); // mol / (cm^3 / s)
        }
    }

    // Add contributions from each reaction for total rate w_i
    for (int i = 0; i < RS.n_species; ++i) {
        rates[i] = 0.0;
        for (auto& r : RS.reactions) {
            rates[i] += RS.MWs[i] * r.dXdt[i] * 1e6;
        }
    }
}

#endif