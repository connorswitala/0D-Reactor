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


void find_Ts(ReactionSet& RS,
             double E_total,        // J/m^3
             double E_vib,          // J/m^3  (mixture vibrational energy density)
             const double* rho_s,   // kg/m^3
             double* Ts)            // Ts[0]=Ttr, Ts[1]=Tv, Ts[2]=Te (optional)
{
    // ---- density ----
    double rho_tot = 0.0;
    for (int i = 0; i < RS.n_species; ++i) rho_tot += rho_s[i];

    // ---- solve Tv from Ev(Tv) = E_vib ----
    double Tv = 300.0; // initial guess (use previous step if you have it)
    const double Tv_min = 1.0;
    const double Tv_max = 5.0e4;

    for (int iter = 0; iter < 50; ++iter) {
        double Ev_model = 0.0;   // J/m^3
        double dEv_dT   = 0.0;   // J/(m^3 K)

        for (int i = 0; i < RS.n_species; ++i) {
            if (!RS.species[i].mol) continue; // atoms: no vibrational energy in this model

            const double theta = RS.species[i].theta_v;   // K
            const double Ri    = RS.species[i].R;         // J/(kg K)
            const double rhoi  = rho_s[i];

            // guard
            double T = std::max(Tv, Tv_min);
            double x = theta / T;

            // exp(x) can overflow if x is huge; clamp x
            if (x > 700.0) x = 700.0;
            double ex = std::exp(x);
            double denom = ex - 1.0;

            // e_v,i(T) = R_i * theta / (exp(theta/T)-1)
            double ev = Ri * theta / denom;               // J/kg
            Ev_model += rhoi * ev;                        // J/m^3

            // derivative: de_v/dT = R_i * theta^2 / T^2 * exp(x) / (exp(x)-1)^2
            double dev_dT = Ri * (theta * theta) / (T * T) * ex / (denom * denom);  // J/(kg K)
            dEv_dT += rhoi * dev_dT;                      // J/(m^3 K)
        }

        double F  = Ev_model - E_vib;
        double Fp = dEv_dT;

        // convergence on energy residual
        if (std::abs(F) < 1e-10 * std::max(1.0, std::abs(E_vib))) break;

        // avoid divide-by-zero
        if (Fp <= 0.0) break;

        double dT = -F / Fp;
        // optional damping if needed
        Tv += dT;

        Tv = std::min(std::max(Tv, Tv_min), Tv_max);

        if (std::abs(dT) < 1e-8 * Tv) break;
    }

    Ts[1] = Tv;

    // ---- compute Ttr from Etr = E_total - E_vib - sum rho_i ef_i ----
    double Etr = E_total - E_vib;  // J/m^3

    for (int i = 0; i < RS.n_species; ++i)
        Etr -= rho_s[i] * RS.species[i].ef;   // ef [J/kg] consistent with E_total definition

    double sum_rho_cv = 0.0; // J/(m^3 K)
    for (int i = 0; i < RS.n_species; ++i) {
        double rot = RS.species[i].mol ? 1.0 : 0.0; // diatomic: +1, monatomic: +0
        sum_rho_cv += rho_s[i] * (1.5 + rot) * RS.species[i].R;
    }

    double Ttr = Etr / sum_rho_cv;
    Ts[0] = Ttr;

    // Ts[2] if you have Te model; otherwise leave as is
}


void landau_teller(ReactionSet& RS, 
                   double* Ts, 
                   double* rho_s, 
                   double& p, 
                   double dt, 
                   double& Q) {
    // Placeholder for Landau-Teller vibrational relaxation implementation
    // This function would update Tv based on Landau-Teller equations

    double Tv = Ts[1];
    double Tr = Ts[0];
    p /= 101325.0; // Convert to atm for the formula

    std::vector<double> X(RS.n_species);
    mass_to_mole(RS, rho_s, X.data());

    std::vector<double> tau_sr(RS.n_species * RS.n_species, 0.0); 
    std::vector<double> tau_v(RS.n_species, 0.0);
    std::vector<double> A(RS.n_species * RS.n_species, 0.0);
    std::vector<double> B(RS.n_species * RS.n_species, 0.0);


    double a = pow(10.0, -4.5);
    a *= 1.16;

    double b = pow(10.0, -0.75);
    b *= 0.015;

    double c;

    for (int i = 0; i < RS.n_species; ++i) {
        for (int j = 0; j < RS.n_species; ++j) {

            if (RS.species[i].mol == false) 
                continue;

            int idx = i * RS.n_species + j;
            c = (RS.species[i].mw * RS.species[j].mw) / (RS.species[i].mw + RS.species[j].mw);
            A[idx] = a * sqrt(c) * RS.species[i].theta_v;       // Placeholder value
            B[idx] = b * pow(c, 0.25);                          // Placeholder value
        }
    }

    for (int i = 0; i < RS.n_species; ++i) {
        for (int j = 0; j < RS.n_species; ++j) {

            if (RS.species[i].mol == false) 
                continue;

            int idx = i * RS.n_species + j;
            tau_sr[idx] = 1 / p * exp(A[idx] * (pow(Ts[0], -1.0/3.0) - B[idx])- 18.42);
        }
    }

    for (int i = 0; i < RS.n_species; ++i) {
        if (RS.species[i].mol == false) 
            continue;

        double deno = 0.0;
        double num = 0.0;

        for (int j = 0; j < RS.n_species; ++j) {

            if (RS.species[j].mol == false) 
                continue;

            int idx = i * RS.n_species + j;

            deno += X[j] / tau_sr[idx];
            num += X[j];
        }   
        tau_v[i] = num / deno;
    }

    double Qs;
    for (int i = 0; i < RS.n_species; ++i) {
        if (RS.species[i].mol == false) 
            continue;

        double ev = RS.species[i].R * RS.species[i].theta_v / (exp(RS.species[i].theta_v / Ts[1]) - 1.0);
        double ev_r = RS.species[i].R * RS.species[i].theta_v / (exp(RS.species[i].theta_v / Ts[0]) - 1.0);

        Qs = rho_s[i] * (ev_r - ev) / tau_v[i];
        Q += Qs;
    }
}


#endif