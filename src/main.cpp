#include "reader.h"
#include "chemistry.h"
#include "energy.h"


inline void write_0D_header(
    const std::string& filename,
    const ReactionSet RS
    ) {

    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out) throw std::runtime_error("Could not open " + filename);

    out << "t,T_tr,T_v";
    for (const auto& nm : RS.species) out << ",rho_" << nm.name;
    out << "\n";
}   



inline void append_0D_row(
    const std::string& filename,
    double t,
    const double* Ts,        // Ts[0]=T_tr, Ts[1]=T_v
    const double* rho_s,
    int n_species
    ) {
        std::ofstream out(filename, std::ios::out | std::ios::app);
        if (!out) throw std::runtime_error("Could not open " + filename + " for append");

        out << t << "," << Ts[0] << "," << Ts[1];
        for (int i = 0; i < n_species; ++i) out << "," << rho_s[i];
        out << "\n";
    }

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/park-90-rates.xml";
    std::string species_data = "../data/species-data.xml";

    read_species_data(species_data, rxns);
    read_rates(park90, rxns);

    std::vector<double> rates(rxns.n_species, 0.0);
    std::vector<double> rho_s(rxns.n_species, 0.0);

    double Ts[3] = {15000.0, 300.0, 1000.0};   // T_tr, T_v, T_e

    double P;
    double E;
    double Ev;
    double rho = 0.0;
    double R_mix;
    double Q; 

    std::string filename = "../files/0DReactor.csv";

    // Initialize species densities (kg/m^3)
    for (int i = 0; i < rxns.n_species; ++i) {
        if (rxns.species[i].name == "O2") 
            rho_s[i] = 0.233; // Example density for O2
        else if (rxns.species[i].name == "N2") 
            rho_s[i] = 0.767; // Example density for N2
        else
        rho_s[i] = 0.0; // Arbitrary small density
    }

    mixture_R(rxns, rho_s.data(), R_mix);
    total_density(rxns, rho_s.data(), rho);
    compute_pressure(P, rho, R_mix, Ts[0]);

    // Precompute species specific gas constants and energies of formation
    for (int i = 0; i < rxns.n_species; ++i) {
        rxns.species[i].R = ugconn / rxns.species[i].mw;
        rxns.species[i].ef = rxns.species[i].hf / rxns.species[i].mw - rxns.species[i].R * 298.15; // convert to kJ/kg
    }

    initialize_energy(rxns, E, Ev, rho_s.data(), Ts);

    double t = 0.0;
    double dt = 1e-14; // time step in seconds

    int NWRITE = 1000;

    write_0D_header(filename, rxns);

    int counter = 0;

    while (t < 1e-5) {

        mixture_R(rxns, rho_s.data(), R_mix);
        total_density(rxns, rho_s.data(), rho);
        compute_pressure(P, rho, R_mix, Ts[0]);

        std::cout << "-- Time: " << t << " s, T_tr: " << Ts[0] << " K, T_v: " << Ts[1] << " K, P: " << P << " Pa" << std::endl;

        for (int i = 0; i < rxns.n_species; ++i) {
            std::cout << "-- " << rxns.species[i].name << " Density: " << rho_s[i] << " kg/m^3" << std::endl;
        }

        find_Ts(rxns, E, Ev, rho_s.data(), Ts);
        compute_rates(rates.data(), rxns, rho_s.data(), Ts, P);        
        landau_teller(rxns, Ts, rho_s.data(), P, dt, Q);

        for (int i = 0; i < rxns.n_species; ++i) {
            rho_s[i] += rates[i] * dt;
            if (rho_s[i] < 0.0) 
                std::cout << "Negative density found!" << rxns.species[i].name << rho_s[i] << std::endl; // Prevent negative densities
            Ev += Q * dt;     
        }
        t += dt;     
        
        
        if (counter % NWRITE == 0) {
            append_0D_row(filename, t, Ts, rho_s.data(), rxns.n_species);
        }

    }

    std::cout << "DONE" << std::endl;


    return 0;

}
