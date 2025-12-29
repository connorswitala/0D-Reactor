#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "../xml/reader.h"
#include "dataStructures.h"

int find_nd_level(const std::vector<double>& Ns, double n) {
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

void compute_rates(double* rates, // Production rates array
    ReactionSet& RS, // Reaction set
    double* rho_s,  // Species densities
    double* Ts,     // Temperatures
    double& P)      // Pressure
    {

    // Molar concentration array
    std::vector<double> C(RS.n_species);

    // Compute concentrations
    for (int i = 0; i < RS.n_species; ++i) 
        C[i] = rho_s[i] / (RS.MWs[i] * 1000.0); // mol/m^3

    for (auto& r : RS.reactions) {

        double T, gamma, n;

        std::vector<double> A(r.Keq_N);
        n = P/(boltzmann * Ts[0] * 1e6);    // Compute number density for curve fits !NEEDS CHANGE
        int level = find_nd_level(r.Ns, n); // Find number density level
 
        // Find Keq coefficients by linear interpolation
        for (int i = 0; i < r.Keq_N; ++i) {
            A[i] = lerp(r.Keq_As[level * r.Keq_N + i], 
                r.Ns[level], 
                r.Keq_As[(level + 1) * r.Keq_N + i], 
                r.Ns[level + 1], 
                n);
        }

        // Define rate-controlling temperature
        if (r.temp == "Ta") {
            T = sqrt(Ts[0] * Ts[1]);
        }
        else if (r.temp == "Te") {
            T = Ts[2];
        }
        else {
            T = Ts[0];
        }

        double exp_term = A[0] * (T / 10000.0)
                            + A[1]
                            + A[2] * log(10000.0 / T)
                            + A[3] * (10000.0 / T)
                            + A[4] * (10000.0 * 10000.0) / (T * T);

        double Keq = exp(exp_term); // Equilibrium constant

        double R, f = 1.0, b = 1.0, nu;
        gamma = r.C * pow(T, r.N) * exp(-r.Ea/T);
            
        for (int i = 0; i < RS.n_species; ++i) {
            nu = RS.nus_f[RS.n_species * r.id + i];
            if (nu == 0)
                continue;
            f *= pow(C[i], nu);
        }

        for (int i = 0; i < RS.n_species; ++i) {
            nu = RS.nus_b[RS.n_species * r.id + i];
            if (nu == 0)
                continue;
            b *= pow(C[i], nu);
        }

        if (r.third_body) {
            double M_eff = 0.0;
            for (int i = 0; i < RS.n_species; ++i) {
                M_eff += C[i] * r.efficiencies[i];
            }

            R = gamma * M_eff * (f - b / Keq);
        }
        else {
            R = gamma * (f - b / Keq);
        }    

        for (int i = 0; i < RS.n_species; ++i) {
            int nu_f = RS.nus_f[RS.n_species * r.id + i];
            int nu_b = RS.nus_b[RS.n_species * r.id + i];

            r.dXdt[i] = R * (nu_b - nu_f);
        }
    }

    for (int i = 0; i < RS.n_species; ++i) {
        rates[i] = 0.0;
        for (auto& r : RS.reactions) {
            rates[i] += RS.MWs[i] * 1000 * r.dXdt[i];
        }
    }
}

#endif