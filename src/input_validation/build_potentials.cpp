#include "input_validation/validate.h"
#include "utility/logger.h"

#include "potentials/coulomb_pot.h"
#include "potentials/yukawa_pot.h"
#include "potentials/exponential_pot.h"

Potential::Ptr_t BuildPotential(const nlohmann::json& potential_item) {
    // make potential object
    if (ToLower(potential_item["type"]) == "coulomb") {
        double Z = potential_item["Z"];
        double x = 0, y = 0, z = 0;

        if (potential_item.contains("location")) {
            x = potential_item["location"][0];
            y = potential_item["location"][1];
            z = potential_item["location"][2];
        }
        
        CoulombPotential* coulomb = new CoulombPotential(Z, x, y, z);
        return Potential::Ptr_t(coulomb);
    } else if (ToLower(potential_item["type"]) == "yukawa") {
        double Z = potential_item["Z"];
        double D = potential_item["decay"];
        double x = 0, y = 0, z = 0;

        if (potential_item.contains("location")) {
            x = potential_item["location"][0];
            y = potential_item["location"][1];
            z = potential_item["location"][2];
        }
        YukawaPotential* yukawa = new YukawaPotential(Z, D, x, y, z);
        return Potential::Ptr_t(yukawa);
    } else if (ToLower(potential_item["type"]) == "exponential") {
        double Z = potential_item["amplitude"];
        double D = potential_item["decay"];
        double x = 0, y = 0, z = 0;

        if (potential_item.contains("location")) {
            x = potential_item["location"][0];
            y = potential_item["location"][1];
            z = potential_item["location"][2];
        }
        ExponentialPotential* expon = new ExponentialPotential(Z, D, x, y, z);
        return Potential::Ptr_t(expon);
    }
    Log::debug(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    return nullptr;
}