#include "input_validation/validate.h"


bool ValidatePotentials(const nlohmann::json& input) {
    if (!input.contains("potentials")) return true;            // no need to check anything since potentials are optional
    if (input.contains("potentials") && !input["potentials"].is_array()) {
        Log::critical("Optional entry \"potentials\" must be an array.");
        return false;
    }

    auto& potentials = input["potentials"];
    for (auto& term : potentials) {
        if (!(term.contains("type") && term["type"].is_string())) {
            MustContain("type", "string");
            return false;
        }

        if (ToLower(term["type"]) == "coulomb") {
            if (!(term.contains("Z") && term["Z"].is_number())) {
                MustContain("Z", "number", "Coulomb term");
                return false;
            }
            if (term.contains("location") && !term["location"].is_array()) {
                Log::critical("Optional entry \"location\" must be a vector in Coloumb term.");
                return false;
            }
        } else if (ToLower(term["type"]) == "yukawa") {
            if (!(term.contains("Z") && term["Z"].is_number())) {
                MustContain("Z", "number", "Yukawa term");
                return false;
            }
            if (!(term.contains("decay") && term["decay"].is_number())) {
                MustContain("decay", "(positive) number", "Yukawa term");
                return false;
            }
            if (term.contains("location") && !term["location"].is_array()) {
                Log::critical("Optional entry \"location\" must be a vector in Yukawa term.");
                return false;
            }
        } else if (ToLower(term["type"]) == "exponential") {
            if (!(term.contains("amplitude") && term["amplitude"].is_number())) {
                MustContain("amplitude", "number", "exponential term");
                return false;
            }
            if (!(term.contains("decay") && term["decay"].is_number())) {
                MustContain("decay", "(positive) number", "exponential term");
                return false;
            }
        }
    }

    return true;
}