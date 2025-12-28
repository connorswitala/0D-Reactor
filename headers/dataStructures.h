#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include "common.h"

struct Species {
    std::string name;           // Name of species
    std::vector<double> lewis;  // NASA polynomial coefficients
    int q;                      // Charge of species
    double hf, ef;              // Enthalpy of formation
    double href;                // Reference enthalpy 
    bool mol;                   // Is molecule boolean
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

    std::vector<int> sp_in_rxn;            // Flags for species

    std::vector<std::string> reactants;     // Reactants in reaction
    std::vector<std::string> products;      // Products in reaction
    std::vector<double> efficiencies;       // Efficiencies for 3rd body
    std::vector<double> Keq_As;             // Equilibrium constant coefficients
    
    int Keq_N;                              // Number of Keq coefficients

    std::vector<double> Ns;                 // Number density levels    
    int nd_levels;                          // Number of density levels

    std::vector<double> dXdt;               // Rate of change of species concentrations
    double Er;                              // Energy of reaction
    
};

struct ReactionSet{
    std::vector<Species> species;           // List of species in reaction set
    std::vector<Reaction> reactions;        // List of reactions
    std::vector<int> nus_f;                 // Total forward reaction stoichiometric set
    std::vector<int> nus_b;                 // Total backward reaction stoichiometric setIDs  
    std::vector<double> MWs;                // Molecular weights of species
    std::vector<double> Rs;                 // Specific gas constants of species
    std::vector<double> theta_vs;            // Characteristic vibrational temperatures of species
    int n_reactions;                        // Number of reactions
    int n_species;                          // Number of species

};


#endif