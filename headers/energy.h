#ifndef ENERGY_H
#define ENERGY_H

#include "dataStructures.h"
#include "common.h"

void initialize_energy(double& E_tot, 
                       double& E_vib, 
                       const ReactionSet& RS,
                       const double* rho_s, 
                       const double* Ts) {

    E_tot = 0.0;
    E_vib = 0.0;

    for (int i = 0; i < RS.n_species; ++i) {

        // Add translational energy
        E_tot += rho_s[i] * 3.0 / 2.0 * RS.Rs[i] * Ts[0];

        // Add rotational energy if molecule
        if (RS.species[i].mol) {
            E_tot += rho_s[i] * RS.Rs[i] * Ts[0];
        }

        // Add vibrational energy if molecule
        if (RS.species[i].mol) {
            double theta_v = RS.theta_vs[i];
            double ev = RS.Rs[i] * theta_v / (exp(theta_v / Ts[1]) - 1.0);
            E_tot += rho_s[i] * ev;
            E_vib += rho_s[i] * ev;
        }

        // Add formation energy 
        E_tot += rho_s[i] * RS.species[i].ef;
    }
}


void find_Ts(double* Ts,
             const ReactionSet& RS,
             const double E_total,        // J/m^3
             const double E_vib,          // J/m^3  (mixture vibrational energy density)
             const double* rho_s,   // kg/m^3
             const double& rho_tot) {          // Ts[0]=Ttr, Ts[1]=Tv, Ts[2]=Te (optional)

    // ---- solve Tv from Ev(Tv) = E_vib ----
    double Tv = Ts[1]; // initial guess (use previous step if you have it)

    for (int iter = 0; iter < 50; ++iter) {
        double Ev_model = 0.0;   // J/m^3
        double dEv_dT   = 0.0;   // J/(m^3 K)

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
            Ev_model += rhoi * ev;                        // J/m^3

            // derivative: de_v/dT = R_i * theta^2 / T^2 * exp(x) / (exp(x)-1)^2
            double dev_dT = Ri * (theta * theta) / (Tv * Tv) * ex / (denom * denom);  // J/(kg K)
            dEv_dT += rhoi * dev_dT;                      // J/(m^3 K)
        }

        double F  = Ev_model - E_vib;
        double Fp = dEv_dT;

        // convergence on energy residual
        if (std::abs(F) < 1e-10 * std::max(1.0, std::abs(E_vib))) 
            break;

        // avoid divide-by-zero
        if (Fp <= 0.0) 
            break;

        double dT = -F / Fp;
        // optional damping if needed
        Tv += dT;

        if (std::abs(dT) < 1e-8 * Tv) 
            break;
    }

    Ts[1] = Tv;

    // ---- compute Ttr from Etr = E_total - E_vib - sum rho_i ef_i ----
    double Etr = E_total - E_vib;  // J/m^3

    for (int i = 0; i < RS.n_species; ++i)
        Etr -= rho_s[i] * RS.species[i].ef;   // ef [J/kg] consistent with E_total definition

    double sum_rho_cv = 0.0; // J/(m^3 K)
    for (int i = 0; i < RS.n_species; ++i) {
        double rot = RS.species[i].mol ? 1.0 : 0.0; // diatomic: +1, monatomic: +0
        sum_rho_cv += rho_s[i] * (1.5 + rot) * RS.Rs[i];
    }

    double Ttr = Etr / sum_rho_cv;
    Ts[0] = Ttr;
}

void landau_teller(double& Q,
                   ReactionSet& RS,
                   const double* Ts,
                   const double* rho_s,
                   const std::vector<double>& Ys,
                   const std::vector<double>& Xs,
                   const double p_Pa) {

    Q = 0.0;

    const double T  = Ts[0];
    const double Tv = Ts[1];
    const double p_atm = p_Pa / 101325.0;
    
    if (p_atm <= 0.0) 
        return;

    std::vector<double> X(RS.n_species);
    mass_to_mole_frac(X.data(), Ys.data(), RS.MWs.data(), RS.n_species);

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

            const double mws = RS.MWs[s]; // confirm units!
            const double mwr = RS.MWs[r];

            const double c = (mws * mwr) / (mws + mwr); // reduced "mass"

            const double A = a * std::sqrt(c) * std::pow(RS.theta_vs[s], 4.0/3.0);
            const double B = b * std::pow(c, 0.25);

            // MW-like form (as you coded)
            tau_sr[idx] = (1.0 / p_atm) * std::exp(A * (std::pow(T, -1.0/3.0) - B) - 18.42);

            n = p_Pa * X[r] / (boltzmann * T); // number density of collider j
            tau_sr[idx] += 1.0 / (C * sigma * n);   // Parks correction
        }
    }

    // Mixture-averaged tau_v for each molecular species i
    for (int s = 0; s < RS.n_species; ++s) {

        if (!RS.species[s].mol) 
            continue;

        double deno = 0.0;
        double num  = 0.0;

        for (int r = 0; r < RS.n_species; ++r) {
            
            if (!RS.species[r].mol) 
                continue;            

            int idx = s * RS.n_species + r;            

            deno += X[r] / tau_sr[idx];
            num  += X[r];
        }

        tau_v[s] = num / deno;
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


#endif