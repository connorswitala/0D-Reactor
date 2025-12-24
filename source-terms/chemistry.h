#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "common.h"
#include "XMLreader.h"

struct Species {
    std::string name;
    double mw;
};

struct Reaction{
    double kf;                              // Forward rate coefficient
    double Keq;                             // Equilibrium constant
    std::vector<double> Keq_coeff;          // Coefficients for K_eq
    bool thirdBodyReaction;                 // Third-body reaction boolean
    std::vector<double> efficiencies;       // Efficiencies for 3rd body
    double eta;                             // Temperature
    double C;                               // Leading coefficient
    double Ea;                              // Activation energy
    std::string index;
    std::vector<std::string> reactants;
    std::vector<std::string> products;
    std::string equation;
    bool ionized;
};

struct ReactionSet{
    std::vector<Species> species;
    std::vector<Reaction> reactions;
    std::vector<int> nus_f;
    std::vector<int> nus_b;
    int n_reactions;
    int n_species;
};


ReactionSet read_rates(std::string& filename) {

    ReactionSet RS;

    try {
        std::string xml = read_file(filename);
        xmlNode root = parse_document(xml);

        // Find root name
        if (root.name != "chemistry")
            throw std::runtime_error("-- Expected <chemistry> root");

        std::string version = require_attr(root, "version");
        std::cout << "-- Chemistry version = " << version << "\n";
        

        // Find species set
        const xmlNode& ss = require_child(root, "species_set");

        for (const auto& ch : ss.children) {
       
            if (ch.name != "species") 
                continue;

            Species sp;
            sp.name = require_attr(ch, "name");
            sp.mw = to_double(require_attr(ch, "MW")); 

            RS.species.push_back(std::move(sp));
        }

        RS.n_species = RS.species.size();

        std::cout << "-- Species count = " << RS.n_species << "\n";
        std::cout << "-- Molecular weights: " << "\n";

        // Print species and molecular weights
        for (auto& s : RS.species) 
            std::cout << "\t" << s.name << " = " << s.mw << "\n";

        // Enter <reactions>
        const xmlNode& rxns = require_child(root, "reactions");

        // Read how many reactions there are
        if (const xmlNode* num = rxns.child("number")) {
            RS.n_reactions = to_int(num->text);
            std::cout << "-- Number of reactions = " << RS.n_reactions << "\n";
        } else {
            throw std::runtime_error("Missing <number> under reactions");
        }

        // // Allocate sized for stoichiometric coefficient arrays
        // RS.nus_f = std::vector<int>(RS.n_reactions * RS.n_species, 0);
        // RS.nus_b = std::vector<int>(RS.n_reactions * RS.n_species, 0);
        

        // // Go through <reaction> blocks
        // for (const auto& rnode : rxns.children) {

        //     if (rnode.name != "reaction") 
        //         continue;

        //     Reaction R;
        //     int id;
        //     id = to_int(require_attr(rnode, "id"));
        //     R.ionized = to_bool(require_attr(rnode, "ionized"));
        //     R.equation = require_child(rnode, "equation").text; 

        //     const xmlNode& react = require_child(rnode, "reactants");

        //     for (auto& sp : react.children) {
        //         if (sp.name != "sp")
        //             continue;

        //         std::string name;
        //         double nu;

        //         name = require_attr(sp, "name");
        //         nu = to_double(require_attr(sp, "nu"));

                
        //         for (int i = 0; i < RS.n_species; ++i)
                

        //     }


        // }
        


    } catch (const std::exception& e) {
        std::cerr << "Parse error: " << e.what() << "\n";
    }

};


#endif