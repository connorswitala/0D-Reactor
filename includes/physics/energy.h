#ifndef ENERGY_H
#define ENERGY_H

#include "dataStructures.h"
#include "../common/common.h"

/**
 *      Functions:
 * 
 *      Ev_sho() - form vibration energy using simple harmonic oscillator
 * 
 *      invert_Ev() - Find T_v from Ev. Only works with single species
 * 
 *      
 */


void Ev_sho(double& Ev,                     // Vibrational energy [J/kg] (to be calculated)
            const double& Tv,               // Vibrational temperature [K]
            const double& rho_s,            // Density of species [kg/m^3]
            const double& R_s,              // Speciific gas constant [J/kg-K]
            const double& theta_v) {        // Characteristic vibrational temperature [K]

    Ev = rho_s * R_s * theta_v / (exp(theta_v / Tv) - 1.0);
}

void invert_Ev(double& Tv,              // Vibrational temperature [K] (to be calculated)
            double& Ev,                 // Vibrational energy [J/kg]
            double& rho_s,              // Partial densities [kg/m^3]
            double& R_s,                // Speciific gas constant [J/kg-K]
            double& theta_v) {          // Characteristic vibrational temperature [K]

    Tv = theta_v / (log(theta_v * R_s * rho_s / Ev + 1));
}

void initialize_energy(double& E_tot,           // Total energy [J/m^3] (to be calculated)
                       double& E_vib,           // Total vibration energy [J/m^3] (to be calculated)
                       const ReactionSet& RS,   // Reaction set
                       const double* rho_s,     // Partial densities [kg/m^3]
                       const double* Ts) {      // Temperatures [K]

    E_tot = 0.0;
    E_vib = 0.0;
    double ev;

    for (int i = 0; i < RS.n_species; ++i) {

        // Add translational energy
        E_tot += rho_s[i] * 3.0 / 2.0 * RS.Rs[i] * Ts[0];

        // Add rotational energy if molecule
        if (RS.species[i].mol) {
            E_tot += rho_s[i] * RS.Rs[i] * Ts[0];
        }

        // Add vibrational energy if molecule
        if (RS.species[i].mol) {
            Ev_sho(ev, Ts[1], rho_s[i], RS.Rs[i], RS.theta_vs[i]);
            E_tot += ev;
            E_vib += ev;
        }

        // Add formation energy 
        E_tot += rho_s[i] * RS.efs[i];
    }
}

void find_Ts(double* Ts,                // Temperatures [K] (to be calculated)
             const ReactionSet& RS,     // Reaction Set
             const double E_total,      // Total energy [J/m^3]
             const double E_vib,        // Vibrational energy [J/m^3]
             const double* rho_s,       // Partial densities [kg/m^3]
             const double& rho_tot) {   // Totaal density [kg/m^3]

    double Tv = Ts[1]; // initial guess

    // Newton solver for Tv
    for (int iter = 0; iter < 50; ++iter) {

        double F = 0.0;   // function = 0
        double Fp = 0.0;  // derivative of function

        for (int i = 0; i < RS.n_species; ++i) {
            if (!RS.species[i].mol) 
                continue; // atoms: no vibrational energy in this model

            const double theta = RS.theta_vs[i];   // K
            const double Ri    = RS.Rs[i];         // J/(kg K)
            const double rhoi  = rho_s[i];

            double x = theta / Tv;

            // exp(x) can overflow if x is huge; clamp x
            if (x > 700.0) x = 700.0;
            double ex = std::exp(x);
            double denom = ex - 1.0;

            // e_v,i(T) = R_i * theta / (exp(theta/T)-1)
            double ev = Ri * theta / denom;               // J/kg
            F += rhoi * ev;                        // J/m^3

            // derivative: de_v/dT = R_i * theta^2 / T^2 * exp(x) / (exp(x)-1)^2
            double dev_dT = Ri * (theta * theta) / (Tv * Tv) * ex / (denom * denom);  // J/(kg K)
            Fp += rhoi * dev_dT;                      // J/(m^3 K)
        }

        F -= E_vib;

        // convergence on energy residual
        if (std::abs(F) < 1e-10) 
            break;

        Tv -= F / Fp;
    }

    Ts[1] = Tv;

    // ---- compute Ttr from Etr = E_total - E_vib - sum rho_i ef_i ----
    double Etr = E_total - E_vib;  // J/m^3

    for (int i = 0; i < RS.n_species; ++i)
        Etr -= rho_s[i] * RS.efs[i];   // ef [J/kg]

    double sum_rho_cv = 0.0; // J/(m^3 K)
    for (int i = 0; i < RS.n_species; ++i) {
        double rot = RS.species[i].mol ? 1.0 : 0.0; // diatomic: +1, monatomic: +0
        sum_rho_cv += rho_s[i] * (1.5 + rot) * RS.Rs[i];
    }

    double Ttr = Etr / sum_rho_cv;
    Ts[0] = Ttr;
}

void landau_teller(double& Q,                       // Source term (to be calculated)
                   ReactionSet& RS,                 // Reaction set
                   const double* Ts,                // Temperatures [K]
                   const double* rho_s,             // Partial densities [kg/m^3]
                   const std::vector<double>& Ys,   // Mass fractions
                   const std::vector<double>& Xs,   // Molar fractions
                   const double p_Pa) {             // Pressure [Pa]

    Q = 0.0;

    const double T  = Ts[0];    // T_translational-rotational
    const double Tv = Ts[1];    // T_vibrational-electronic
    const double p_atm = p_Pa / 101325.0; // Pressure [atm]

    std::vector<double> tau_sr(RS.n_species * RS.n_species, 0.0);
    std::vector<double> tau_v (RS.n_species, 0.0);

    // constants: make sure these match the correlation + units you want
    double a = std::pow(10.0, -4.5) * 1.16;
    double b = std::pow(10.0, -0.75) * 0.015;
    double C, sigma, n;

    // Build pair relaxation times for molecular relaxers i with ANY collider j
    for (int s = 0; s < RS.n_species; ++s) {

        if (!RS.species[s].mol) 
            continue;

        C = std::sqrt(8.0 * RS.Rs[s] * T / pi); // mean thermal speed
        sigma = 3e-21 * 50000.0 * 50000.0 / (T * T); // collision cross section

        for (int r = 0; r < RS.n_species; ++r) {

            int idx = s * RS.n_species + r;

            const double mws = RS.MWs[s] * 1000.0; // confirm units!
            const double mwr = RS.MWs[r] * 1000.0;

            const double c = (mws * mwr) / (mws + mwr); // reduced "mass"

            const double A = a * std::sqrt(c) * std::pow(RS.theta_vs[s], 4.0/3.0);
            const double B = b * std::pow(c, 0.25);

            // MW-like form (as you coded)
            double tau_mw = (1.0 / p_atm) * std::exp(A * (std::pow(T, -1.0/3.0) - B) - 18.42);

            n = p_Pa * Xs[r] / (boltzmann * T); // number density of collider j
            double tau_p = 1.0 / (C * sigma * n);   // Parks correction

            tau_sr[idx] = tau_mw + tau_p;

        }
    }

    // Mixture-averaged tau_v for each molecular species i
    for (int s = 0; s < RS.n_species; ++s) {

        if (!RS.species[s].mol) 
            continue;

        double sum = 0.0;

        for (int r = 0; r < RS.n_species; ++r) {           

            int idx = s * RS.n_species + r;
            sum += Xs[r] / tau_sr[idx];
        }

        tau_v[s] = 1.0 / sum;
    }

    // Q_VT = sum rho_i (e_v*(T) - e_v(Tv)) / tau_v,i
    for (int i = 0; i < RS.n_species; ++i) {

        if (!RS.species[i].mol) 
            continue;

        auto ev_sho = [&](double TT){
            double x = RS.theta_vs[i] / TT;
            return RS.Rs[i] * RS.theta_vs[i] / (std::exp(x) - 1.0); // J/kg
        };

        const double ev   = ev_sho(Tv);
        const double eveq = ev_sho(T);

        Q += rho_s[i] * (eveq - ev) / tau_v[i];
    }
}

void chemical_vibrational(double& Q,
                          ReactionSet& RS,
                          double& Tv,
                          std::vector<double> rates) {

    Q = 0.0;

    for (int i = 0; i < RS.n_reactions; ++i) {
        if (!RS.species[i].mol) 
            continue;

        auto ev_sho = [&](double TT){
            double x = RS.theta_vs[i] / TT;
            return RS.Rs[i] * RS.theta_vs[i] / (std::exp(x) - 1.0); // J/kg
        };

        Q += rates[i] * ev_sho(Tv);
    }

}

#endif