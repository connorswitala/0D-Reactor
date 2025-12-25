#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "common.h"
#include "XMLreader.h"

struct Species {
    std::string name;
    double mw;
};

struct Reaction{
    std::vector<double> Keq_coeff;          // Coefficients for K_eq
    bool third_body;                 // Third-body reaction boolean
    std::vector<double> efficiencies;       // Efficiencies for 3rd body
    double eta;                             // Temperature
    double C;                               // Leading coefficient
    double Ea;                              // Activation energy
    double N;                               // Temperature exponent
    std::string temp;                       // Reaction temperature
    int id;                                 // Reaction ID
    std::vector<std::string> reactants;     // Reactants in reaction
    std::vector<std::string> products;      // Products in reaction
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

        // Allocate sized for stoichiometric coefficient arrays
        RS.nus_f = std::vector<int>(RS.n_reactions * RS.n_species, 0);
        RS.nus_b = std::vector<int>(RS.n_reactions * RS.n_species, 0);
        

        // Go through <reaction> blocks
        for (const auto& rnode : rxns.children) {

            if (rnode.name != "reaction") 
                continue;

            Reaction R;
            R.id = to_int(require_attr(rnode, "id"));
            R.ionized = to_bool(require_attr(rnode, "ionized"));
            R.equation = require_child(rnode, "equation").text; 

            const xmlNode& reactants = require_child(rnode, "reactants");

            for (auto& sp : reactants.children) {
                if (sp.name != "sp")
                    continue;

                std::string name;
                double nu;

                name = require_attr(sp, "name");
                nu = to_int(require_attr(sp, "nu"));
                
                for (int i = 0; i < RS.n_species; ++i) {
                    if (RS.species[i].name == name) 
                        RS.nus_f[RS.n_species * R.id + i] = nu;
                }
            }

            const xmlNode& products = require_child(rnode, "products");

            for (auto& sp : products.children) {
                if (sp.name != "sp")
                    continue;

                std::string name;
                int nu;

                name = require_attr(sp, "name");
                nu = to_int(require_attr(sp, "nu"));
                
                for (int i = 0; i < RS.n_species; ++i) {
                    if (RS.species[i].name == name) {
                        RS.nus_b[RS.n_species * R.id + i] = nu;
                        break;
                    }
                }
            }

            const xmlNode& rate = require_child(rnode, "rate");

            // <C units="...">2.0e15</C>
            const xmlNode& Cnode = require_child(rate, "C");

            std::string C_units = require_attr(Cnode, "units");
            if (C_units != "m^3/(mol-sec)")
                throw std::runtime_error("Invalid units for C: " + C_units);

            R.C  = to_double(Cnode.text);

            // <T>Ta</T>
            R.temp = require_child(rate, "T").text;   // string "Ta"

            // <N>-1.5</N>
            R.N  = to_double(require_child(rate, "N").text);

            // <Ea>59500.0</Ea>
            R.Ea = to_double(require_child(rate, "Ea").text);

            if (const xmlNode* third = rnode.child("third-body")) {
                // third-body exists → parse it

                for (const auto& e : third->children) {
                    if (e.name != "eff") 
                        continue;

                    std::string sp = require_attr(e, "sp");
                    double val = to_double(require_attr(e, "efficiency"));
                    
                    for (int i = 0; i < RS.n_species; ++i) {
                        if (RS.species[i].name == sp) {
                            R.efficiencies[i] = val;                            
                        }
                    }
                }

                R.third_body = true;
            } else {
                // no third-body block → bimolecular reaction
                R.third_body = false;
            }

            RS.reactions.push_back(std::move(R));
        }

        for (auto& r : RS.reactions)  {
            std::cout << "Reaction " << r.id << "\t Equation = " << r.equation << ": \n";
            std::cout << "Reactants: ";
            for (int i = 0; i < RS.n_species; ++i) {
                std::cout << RS.species[i].name << " has nu = " << RS.nus_f[RS.n_species * r.id + i] << ", ";
            }
            std::cout << "\nProducts: ";
            for (int i = 0; i < RS.n_species; ++i) {
                std::cout << RS.species[i].name << " has nu = " << RS.nus_b[RS.n_species * r.id + i] << ", ";
            }
            if (r.third_body) {
                std::cout << "\nThird-body efficiencies: ";
                for (int i = 0; i < RS.n_species; ++i) {
                    std::cout << RS.species[i].name << " has efficiency = " << r.efficiencies[i] << ", ";
                }
            }
            std::cout << "\nC = " << r.C << ", T = " << r.temp << ", Ea = " << r.Ea << ", N = " << r.N;
            std::cout << "\n\n";
        }


    } 
    catch (const std::exception& e) {
        std::cerr << "Parse error: " << e.what() << "\n";
    }

    return RS;
}


#endif