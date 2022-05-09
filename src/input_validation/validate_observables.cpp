#include "input_validation/validate.h"

#include <iostream>

bool ValidateObservables(const nlohmann::json& input) {
    if (!input.contains("observables")) return true;            // no need to check anything since observables are optional
    if (input.contains("observables") && !input["observables"].is_object()) {
        Log::critical("Optional entry \"observables\" must be an object.");
        return false;
    }

    auto& observables = input["observables"];
    for (auto& obs_pair : observables.items()) {
        auto& obs = obs_pair.value();
        // might contain 'compute_period'. Defaults to 1.
        if (obs.contains("compute_period") && !obs["compute_period"].is_number()) {
            Log::critical("Optional entry \"compute_period\" must be a number.");
            return false;
        }
        // filename is sometimes optional. If it there make sure it is a number.
        if (obs.contains("filename") && !obs["filename"].is_string()) {
            Log::critical("Optional entry \"filename\" must be a string.");
            return false;
        }
        
        if (obs_pair.key() == "norm") {
            // do not require anything
        } else if (obs_pair.key() == "wavefunction") {
            if (!(obs_pair.value().contains("grid_points") && obs_pair.value()["grid_points"].is_number())) {
                MustContain("grid_points", "number", "wavefunction observable");
                return false;
            }
            // an output file is required
            if (!obs_pair.value().contains("filename")) {
                MustContain("filename", "string", "wavefunction observable");
                return false;
            }
        } else if (obs_pair.key() == "dipole_acc") {
            // an output file is required
            if (!obs_pair.value().contains("filename")) {
                MustContain("filename", "string", "dipole_acc observable");
                return false;
            }
        } else if (obs_pair.key() == "populations") {
            // an output file is required
            if (!obs_pair.value().contains("filename")) {
                MustContain("filename", "string", "populations observable");
                return false;
            }
        }
    }
    return true;
}

