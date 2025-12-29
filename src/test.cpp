#include "../includes/xml/reader.h"

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/N2-park-90-rates.xml";
    std::string species_data = "../data/species-data.xml";

    read_rates(park90, species_data, rxns);

    return 0;
}