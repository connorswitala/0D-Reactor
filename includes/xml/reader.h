#ifndef READER_H
#define READER_H

#include "XMLreader.h"
#include "../physics/dataStructures.h"


static void hr(char ch='-', int w=100) {
    for (int i=0;i<w;i++) std::cout << ch;
    std::cout << "\n";
}

// helper: safe access for nus arrays (prevents accidental OOB in debug)
static int nu_at(const std::vector<int>& nu, int nsp, int rid, int i) {
    const size_t idx = static_cast<size_t>(nsp) * static_cast<size_t>(rid) + static_cast<size_t>(i);
    return nu.at(idx);
}

static std::string trim_one_line(std::string s) {
    // collapse whitespace/newlines for equation printing
    for (auto& c : s) if (c == '\n' || c == '\r' || c == '\t') c = ' ';
    // squeeze multiple spaces
    std::string out;
    out.reserve(s.size());
    bool prev_space = false;
    for (char c : s) {
        bool is_space = (c == ' ');
        if (is_space) {
            if (!prev_space) out.push_back(' ');
        } else {
            out.push_back(c);
        }
        prev_space = is_space;
    }
    // trim ends
    while (!out.empty() && out.front()==' ') out.erase(out.begin());
    while (!out.empty() && out.back() ==' ') out.pop_back();
    return out;
}

void print_reaction_set_pretty(const ReactionSet& RS) {
    constexpr int W = 110;           // overall width
    constexpr int col_sp = 12;       // species name column width
    constexpr int col_val = 10;      // numeric column width

    std::cout << "\n";
    hr('=', W);
    std::cout << "CHEMISTRY SUMMARY\n";
    hr('=', W);
    std::cout << "Species: " << RS.n_species
              << " | Reactions: " << RS.reactions.size()
              << "\n\n";

    // Optional: print species table
    {
        std::cout << "Species (name, MW, theta_v, hf):\n";
        hr('-', W);
        for (int i = 0; i < RS.n_species; ++i) {
            std::cout << "  " << std::left << std::setw(6) << i
                      << std::left << std::setw(col_sp) << RS.species[i].name
                      << "MW= " << std::right << std::setw(10) << std::fixed << std::setprecision(6)
                      << RS.MWs[i]
                      << "  theta_v= " << std::right << std::setw(10) << std::fixed << std::setprecision(4)
                      << RS.theta_vs[i]
                      << "  hf= " << std::right << std::setw(10) << std::fixed << std::setprecision(4)
                      << RS.species[i].hf << "\n";
        }
        hr('-', W);
        std::cout << "\n";
    }

    for (const auto& r : RS.reactions) {
        // Header
        hr('=', W);
        std::cout << "Reaction " << r.id << "\n";
        hr('-', W);

        // Equation (single-line)
        std::cout << "Equation: " << trim_one_line(r.equation) << "\n\n";

        // Rate & equilibrium scalars
        std::cout << "Rate model params:\n";
        std::cout << "  C  = " << std::scientific << std::setprecision(6) << r.C
                  << "   Texp = " << r.Texp
                  << "   Ea = " << std::fixed << std::setprecision(3) << r.Ea
                  << "   N = " << std::setprecision(6) << r.N
                  << "\n";
        std::cout << "  Er = " << std::fixed << std::setprecision(6) << r.Er << "\n\n";

        // Stoichiometry table
        std::cout << "Stoichiometry (nu_f reactants, nu_b products):\n";
        hr('-', W);

        // table header
        std::cout << std::left
                  << std::setw(6) << "idx"
                  << std::setw(col_sp) << "species"
                  << std::right
                  << std::setw(col_val) << "nu_f"
                  << std::setw(col_val) << "nu_b";

        if (r.third_body) {
            std::cout << std::setw(col_val+4) << "eff";
        }
        std::cout << "\n";
        hr('-', W);

        // rows
        for (int i = 0; i < RS.n_species; ++i) {
            const int nu_f = nu_at(RS.nus_f, RS.n_species, r.id, i);
            const int nu_b = nu_at(RS.nus_b, RS.n_species, r.id, i);

            // Only print species that participate OR print all if you prefer.
            if (nu_f == 0 && nu_b == 0 && !r.third_body) continue;

            std::cout << std::left
                      << std::setw(6) << i
                      << std::setw(col_sp) << RS.species[i].name
                      << std::right
                      << std::setw(col_val) << nu_f
                      << std::setw(col_val) << nu_b;

            if (r.third_body) {
                // guard: efficiencies should be sized to n_species
                double eff = (i < (int)r.efficiencies.size()) ? r.efficiencies[i] : 1.0;
                std::cout << std::setw(col_val+4) << std::fixed << std::setprecision(3) << eff;
            }

            std::cout << "\n";
        }

        hr('-', W);
        std::cout << "\n";

        // Equilibrium coefficients (pretty grid)
        std::cout << "Equilibrium constant coefficients (A1..A5 per density level):\n";
        if (r.Keq_As.empty()) {
            std::cout << "  (none)\n";
        } else {
            const int levels = std::max(1, r.nd_levels);  // avoid div-by-zero
            const int cols = 5;                            // A1..A5
            const int rows = (int)r.Keq_As.size() / cols;  // assuming stored as [level][A1..A5]

            // If your storage is different, adjust above.
            std::cout << "  levels = " << levels << ", total coeffs = " << r.Keq_As.size() << "\n\n";

            // Print header
            std::cout << std::left << std::setw(10) << "level"
                      << std::right << std::setw(14) << "A1"
                      << std::setw(14) << "A2"
                      << std::setw(14) << "A3"
                      << std::setw(14) << "A4"
                      << std::setw(14) << "A5"
                      << "\n";
            hr('-', W);

            // Print rows as A1..A5
            for (int lvl = 0; lvl < rows; ++lvl) {
                const int base = lvl * cols;
                if (base + 4 >= (int)r.Keq_As.size()) break;

                std::cout << std::left << std::setw(10) << lvl
                          << std::right << std::setw(14) << std::fixed << std::setprecision(6) << r.Keq_As[base + 0]
                          << std::setw(14) << r.Keq_As[base + 1]
                          << std::setw(14) << r.Keq_As[base + 2]
                          << std::setw(14) << r.Keq_As[base + 3]
                          << std::setw(14) << r.Keq_As[base + 4]
                          << "\n";
            }
        }

        std::cout << "\n";
    }

    hr('=', W);
    std::cout << "End of chemistry summary.\n";
    hr('=', W);
    std::cout << std::flush;
}

void read_species_data(std::string& filename, std::vector<std::string> names, ReactionSet& RS) {

    Species sp;

    try{
        std::string xml = read_file(filename);
        xmlNode root = parse_document(xml);

        // Find root name
        if (root.name != "species_data")
            throw std::runtime_error("-- Expected <species_data> root");

        std::string version = require_attr(root, "version");
        std::cout << "-- Species data version = " << version << "\n";

        int n_species = names.size();
        double mw;
        double theta_v;
        double vibc1;
        bool inlist;
        double Rs;
        double ef;

        // Enter <species> blocks
        for (const auto& ch : root.children) {

            if (ch.name != "species")
                continue;
                
            sp.name = require_child(ch, "name").text;

            inlist = false;
            for (int i = 0; i < names.size(); i++) {
                if (sp.name == names[i]) 
                    inlist = true;
            }

            if (!inlist)
                continue;

            mw = to_double(require_child(ch, "MW").text);
            vibc1 = to_double(require_child(ch, "vibc1").text);

            if (vibc1 > 0.0) {
                theta_v = vibc1 * planck * light_speed / boltzmann * 100.0; // convert to K
                sp.mol = true;
            }
            else {
                theta_v = 0.0;
                sp.mol = false;
            }

            sp.hf = to_double(require_child(ch, "hf").text);

            Rs = ugconn / mw;
            ef = sp.hf / mw - Rs * 298.15;

            // Add to species list
            RS.species.push_back(std::move(sp));
            RS.MWs.push_back(mw);
            RS.theta_vs.push_back(theta_v);
            RS.Rs.push_back(Rs);
            RS.efs.push_back(ef);
            
        }

        RS.n_species = n_species;
        RS.Rs.resize(RS.n_species); 
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Error reading species file '" + filename + "': " + e.what());
    }
}

void read_rates(std::string& rates_file, std::string species_data_file, ReactionSet& RS) {

    try {
        std::string xml = read_file(rates_file);
        xmlNode root = parse_document(xml);

        // Find root name
        if (root.name != "chemistry")
            throw std::runtime_error("-- Expected <chemistry> root");

        std::string version = require_attr(root, "version");
        std::cout << "-- Chemistry version = " << version << "\n";

        // Enter <reactions>
        const xmlNode& rxns = require_child(root, "reactions");

        // Enter <species>
        const xmlNode& species_names = require_child(rxns, "species");
        std::vector<std::string> names;

        for (const auto& spnames : species_names.children) {
            if (spnames.name != "sp")
                continue;

            std::string name = require_attr(spnames, "name");
            names.push_back(name);             
        }   

        // Read in species data for species in this reaction set
        read_species_data(species_data_file, names, RS);

        // Readnumber of reactions
        const xmlNode& num = require_child(rxns, "number");
        RS.n_reactions = to_int(num.text);

        // Allocate sized for stoichiometric coefficient arrays
        RS.nus_f = std::vector<int>(RS.n_reactions * RS.n_species, 0);
        RS.nus_b = std::vector<int>(RS.n_reactions * RS.n_species, 0);

        // Go through <reaction> blocks
        for (const auto& rnode : rxns.children) {

            if (rnode.name != "reaction")
                continue;

            Reaction R;
            R.dXdt = std::vector<double>(RS.n_species);

            // Read in reaction id, ionization, equation
            R.id = to_int(require_attr(rnode, "id"));
            R.ionized = to_bool(require_attr(rnode, "ionized"));
            R.equation = require_child(rnode, "equation").text;

            
            // Read in reaction reactants and stoichiometry      
            const xmlNode& reactants = require_child(rnode, "reactants");

            for (auto& sp : reactants.children) {
                if (sp.name != "sp")
                    continue;

                std::string name;
                double nu;

                name = require_attr(sp, "name");
                nu = to_int(require_attr(sp, "nu"));

                for (int i = 0; i < RS.n_species; ++i) {
                    if (RS.species[i].name == name) {
                        RS.nus_f[RS.n_species * R.id + i] = nu;
                        break;
                    }
                }
            }

            /**
             * Read in reaction products and stoichiometry
             */

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

            /**
             * Read in forward rate coefficients
             */

            const xmlNode& rate = require_child(rnode, "rate");

            const xmlNode& Cnode = require_child(rate, "C");

            std::string C_units = require_attr(Cnode, "units");
            if (C_units != "cm^3/(mol-sec)")
                throw std::runtime_error("Invalid units for C: " + C_units);

            R.C  = to_double(Cnode.text);

            // <T>Ta</T>
            R.Texp = to_double(require_child(rate, "exp").text);

            // <N>-1.5</N>
            R.N  = to_double(require_child(rate, "N").text);

            // <Ea>59500.0</Ea>
            R.Ea = to_double(require_child(rate, "Ea").text);

            /**
             * Read in third body reactions
             */

            R.efficiencies = std::vector<double>(RS.n_species);

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

            /**
             * Read in Energy of reaction
             */

            const xmlNode& eq = require_child(rnode, "equilibrium");

            const xmlNode& Er = require_child(eq, "Er");
            std::string units = require_attr(Er, "units");
            if (units != "J/mol")
                throw std::runtime_error("Incorrect units for energy of reaction: " + units);

            R.Er = to_double(Er.text);

            const xmlNode& levels = require_child(eq, "levels");
            R.nd_levels = to_int(levels.text);

            const xmlNode& ncoeffs = require_child(eq, "ncoeffs");
            R.Keq_N = to_int(ncoeffs.text);

            for (auto& f : eq.children) {
                if (f.name != "fit")
                    continue;

                R.Ns.push_back(to_double(require_attr(f, "N")));
                R.Keq_As.push_back(to_double(require_child(f, "A1").text));
                R.Keq_As.push_back(to_double(require_child(f, "A2").text));
                R.Keq_As.push_back(to_double(require_child(f, "A3").text));
                R.Keq_As.push_back(to_double(require_child(f, "A4").text));
                R.Keq_As.push_back(to_double(require_child(f, "A5").text));


            }
            RS.reactions.push_back(std::move(R));
        }

        print_reaction_set_pretty(RS);
    }
    catch (const std::exception& e) {
        std::cerr << "Parse error: " << e.what() << "\n";
    }
}

#endif