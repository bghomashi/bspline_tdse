#include "input_validation/validate.h"

bool ValidateLasers(const nlohmann::json& input) {
    if (input.contains("lasers") && input["lasers"].is_array()) {
        auto& lasers = input["lasers"];

        for (auto& pulse : lasers) {
            if (!(pulse.contains("cycles_delay") && pulse["cycles_delay"].is_number())) {
                MustContain("cycles_delay", "number");
                return false;
            }
            if (!(pulse.contains("envelope") && pulse["envelope"].is_string())) {
                MustContain("envelope", "string");
                return false;
            }
            if (!(pulse.contains("intensity") && pulse["intensity"].is_number())) {
                MustContain("intensity", "number");
                return false;
            }
            if (!(pulse.contains("num_cycles") && pulse["num_cycles"].is_number())) {
                MustContain("num_cycles", "number");
                return false;
            }
            if (!(pulse.contains("wavelength") && pulse["wavelength"].is_number()) && 
                !(pulse.contains("energy") && pulse["energy"].is_number())) {
                Log::critical("input file must contain wavelength or energy entry.");
                return false;
            }
            if (!(pulse.contains("cep") && pulse["cep"].is_number())) {
                MustContain("cep", "number");
                return false;
            }
            if (pulse.contains("ellipticity") && !pulse["ellipticity"].is_number()) {
                Log::critical("optional parameter \"ellipticity\" must be a number.");
                return false;
            }
            if (!(pulse.contains("polarization_vector") && pulse["polarization_vector"].is_array())) {
                MustContain("polarization_vector", "3-vector", "lasers");
                return false;
            }
            if (!(pulse.contains("poynting_vector") && pulse["poynting_vector"].is_array())) {
                MustContain("poynting_vector", "3-vector", "lasers");
                return false;
            }
        }
    }
    return true;
}