#include "reader.h"
#include "chemistry.h"

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/park-90-rates.xml";

    rxns = read_rates(park90);

    std::vector<double> rates(rxns.n_species, 0.0);
    std::vector<double> rho_s(rxns.n_species, 0.0);
    double Ts[3] = {2500.0, 2500.0, 2500.0};   // T_tr, T_v, T_e
    double P = 345.1623; // Pressure in Pa
    // Initialize species densities (kg/m^3)
    
    for (int i = 0; i < rxns.n_species; ++i) {
        if (rxns.species[i].name == "O2") 
            rho_s[i] = 0.24; // Example density for O2
        else if (rxns.species[i].name == "N2") 
            rho_s[i] = 0.76; // Example density for N2
        else
        rho_s[i] = 0.0; // Arbitrary small density
    }

    compute_rates(rates.data(), rxns, rho_s.data(), Ts, P);

    double sum = 0.0;

    for (int i = 0; i < rxns.n_species; ++i) {
        std::cout << "Rate of production for species " << rxns.species[i].name 
                  << " = " << rates[i] << " kg/m^3-s\n";
        sum += rates[i];
    }

    std::cout << "Sum of rates = " << sum << " kg/m^3-s\n";



    std::cout << "DONE" << std::endl;


    return 0;

}
