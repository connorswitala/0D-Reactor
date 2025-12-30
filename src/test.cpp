#include "../includes/xml/reader.h"
#include "../includes/physics/chemistry.h"
#include "../includes/physics/energy.h"


inline void write_0D_header(
    const std::string& filename,
    const ReactionSet RS
    ) {

    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out) throw std::runtime_error("Could not open " + filename);

    out << "t,T_tr,T_v,rho,p";
    for (const auto& nm : RS.species) out << ",X_" << nm.name;
    out << "\n";
}   

inline void append_0D_row(
    const std::string& filename,
    double t,
    const double* Ts,        // Ts[0]=T_tr, Ts[1]=T_v
    const double* X_s,
    const double& rho,
    const double& p,
    int n_species) {

        std::ofstream out(filename, std::ios::out | std::ios::app);
        if (!out) throw std::runtime_error("Could not open " + filename + " for append");

        out << t << "," << Ts[0] << "," << Ts[1] << "," << rho << "," << p;
        for (int i = 0; i < n_species; ++i) out << "," << X_s[i];
        out << "\n";
    }

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/air5-park-90-rates.xml";
    std::string species_data = "../data/species-data.xml";

    read_rates(park90, species_data, rxns);

    std::vector<double> rates(rxns.n_species, 0.0);
    std::vector<double> rho_s(rxns.n_species, 0.0);
    std::vector<double> Ys(rxns.n_species, 0.0);
    std::vector<double> Xs(rxns.n_species, 0.0);

    double Ts[3] = {15000.0, 1000.0, 1000.0};   // T_tr, T_v, T_e

    double P = 20.42 * 101325.0;
    double E;
    double Ev;
    double rho;
    double R_mix;
    double Q; 

    std::string filename = "../files/0DReactor.csv";

    // Initialize species densities (kg/m^3)
    for (int i = 0; i < rxns.n_species; ++i) {
        if (rxns.species[i].name == "N2") 
            Ys[i] = 0.2; // Example density for O2
        else if (rxns.species[i].name == "O2") 
            Ys[i] = 0.7; // Example density for N2
        else if (rxns.species[i].name == "NO") 
            Ys[i] = 0.05; // Example density for N2
        else if (rxns.species[i].name == "N") 
            Ys[i] = 0.025; // Example density for N2
        else if (rxns.species[i].name == "O") 
            Ys[i] = 0.025; // Example density for N2        
    }

    mixture_gas_constant(R_mix, Ys.data(), rxns.MWs.data(), rxns.n_species);
    rho = P / (R_mix * Ts[0]);
    for (int i = 0; i < rxns.n_species; ++i) {
        rho_s[i] = Ys[i] * rho;
    }
    mass_to_mole_frac(Xs.data(), Ys.data(), rxns.MWs.data(), rxns.n_species);

    // Precompute species specific gas constants and energies of formation
    for (int i = 0; i < rxns.n_species; ++i) {
        rxns.Rs[i] = ugconn / rxns.MWs[i];
        rxns.species[i].ef = rxns.species[i].hf / rxns.MWs[i] - rxns.Rs[i] * 298.15; // convert to J/kg
    }
    
    initialize_energy(E, Ev, rxns, rho_s.data(), Ts);
    find_Ts(Ts, rxns, E, Ev, rho_s.data(), rho);

    std::cout << "Ev = " << Ev << std::endl;
    std::cout << "Tv = " << Ts[1] << std::endl;


    double ev;
    double sum = 0.0;
    for (int i = 0; i < rxns.n_species; ++i) {
        if (!rxns.species[i].mol)
            continue;
        Ev_sho(ev, Ts[1], rho_s[i], rxns.Rs[i], rxns.theta_vs[i]);
        std::cout << rxns.species[i].name << " Ev = " << ev << std::endl;
        sum += ev;
    }

    std::cout << "Sum = " << sum << std::endl;
    std::cout << "DONE" << std::endl;
    return 0;

}
