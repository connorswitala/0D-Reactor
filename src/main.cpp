#include "chemistry.h"

#include <iostream>

int main() {

    ReactionSet rxns;

    std::string park90 = "../data/park-90-rates.xml";

    rxns = read_rates(park90);

    std::cout << "DONE" << std::endl;


    return 0;

}
