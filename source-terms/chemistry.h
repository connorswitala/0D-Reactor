#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "common.h"
#include "XMLreader.h"

struct Species {
    std::string name;                       // Species name
    double mw;                              // Species molecular weight
};

struct Reaction{
    bool third_body;                        // Third-body reaction boolean
    bool ionized;                           // Ionized reaction boolean

    double eta;                             // Temperature
    double C;                               // Leading coefficient
    double Ea;                              // Activation energy
    double N;                               // Temperature exponent
    std::string temp;                       // Reaction temperature

    int id;                                 // Reaction ID
    std::string equation;                   // Equation 

    std::vector<std::string> reactants;     // Reactants in reaction
    std::vector<std::string> products;      // Products in reaction
    std::vector<double> efficiencies;       // Efficiencies for 3rd body
    std::vector<double> Keq_As;             // Equilibrium constant coefficients

    std::vector<double> Ns;
    int nd_levels;
    
    double Er;                              // Energy of reaction
    
};

struct ReactionSet{
    std::vector<Species> species;           // List of species in reaction set
    std::vector<Reaction> reactions;        // List of reactions
    std::vector<int> nus_f;                 // Total forward reaction stoichiometric set
    std::vector<int> nus_b;                 // Total backward reaction stoichiometric set
    int n_reactions;                        // Number of reactions
    int n_species;                          // Number of species
};


void compute_rates(double* rates, 
    ReactionSet RS, 
    double* X, 
    double* Ts, 
    double& P) {

    // Ts[0] = T_tr;
    // Ts[1] = T_v;
    // Ts[2] = T_e;

    for (const auto& r : RS.reactions) {

        double T, gamma, n;

        if (r.temp == "Ta") {

            T = sqrt(Ts[0] * Ts[1]);
            double R = 0.0;

            if (r.third_body) {
                gamma = r.C * pow(T, r.N) * exp(-r.Ea/(bk * T));
                n = P/(bk * Ts[0]);
                
                for (int i; i < RS.n_species; ++i) {
                    


                }



            }
            else {

            }


        } 
        else if (r.temp == "T") {

        }
        else if (r.temp == "Te") {


        }






    }




}





#endif