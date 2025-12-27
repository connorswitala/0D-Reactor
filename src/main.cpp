#include "reader.h"
#include "chemistry.h"
#include "energy.h"

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/park-90-rates.xml";
    std::string species_data = "../data/species-data.xml";

    read_species_data(species_data, rxns);
    read_rates(park90, rxns);

    std::vector<double> rates(rxns.n_species, 0.0);
    std::vector<double> rho_s(rxns.n_species, 0.0);
    double Ts[3] = {300.0, 300.0, 300.0};   // T_tr, T_v, T_e
    double P = 101325.0; // Pressure in Pa
    double E;
    double Ev;
    // Initialize species densities (kg/m^3)

    // Initialize species densities (kg/m^3)
    for (int i = 0; i < rxns.n_species; ++i) {
        if (rxns.species[i].name == "O2") 
            rho_s[i] = 0.24; // Example density for O2
        else if (rxns.species[i].name == "N2") 
            rho_s[i] = 0.76; // Example density for N2
        else
        rho_s[i] = 0.0; // Arbitrary small density
    }

    // Precompute species specific gas constants and energies of formation
    for (int i = 0; i < rxns.n_species; ++i) {
        rxns.species[i].R = ugconn / rxns.species[i].mw;
        rxns.species[i].ef = rxns.species[i].hf / rxns.species[i].mw - rxns.species[i].R * 298.15; // convert to kJ/kg
    }

    initialize_energy(rxns, E, Ev, rho_s.data(), Ts);

    std::cout << "Initial total energy = " << E << " J/kg\n";
    std::cout << "Initial vibrational energy = " << Ev << " J/kg\n";

    compute_rates(rates.data(), rxns, rho_s.data(), Ts, P);

    double sum = 0.0;
    

    // for (int i = 0; i < rxns.n_species; ++i) {
    //     std::cout << "Rate of production for species " << rxns.species[i].name 
    //               << " = " << rates[i] << " kg/m^3-s\n";
    //     sum += rates[i];
    // }

    std::cout << "Sum of rates = " << sum << " kg/m^3-s\n";



    std::cout << "DONE" << std::endl;


    return 0;

}
