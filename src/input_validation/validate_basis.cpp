#include "input_validation/validate.h"

bool ValidateBasis(const nlohmann::json& input) {
    if (!(input.contains("basis") && input["basis"].is_object())) {
        MustContain("basis", "object");
        return false;
    } else {
        auto& basis = input["basis"];
        if (!(basis.contains("order") && basis["order"].is_number())) {
            MustContain("order", "number", "basis");
            return false;
        }
        if (!(basis.contains("node_sequence") && basis["node_sequence"].is_string())) {
            MustContain("node_sequence", "string", "basis");
            return false;
        } 
        if (ToLower(basis["node_sequence"]) == "linear") {
            // fine
        } else if (ToLower(basis["node_sequence"]) == "exponential") {
            if (!(basis.contains("parameter") && basis["parameter"].is_number())) {
                MustContain("parameter", "number", "basis");
                return false;
            }
        } else if (ToLower(basis["node_sequence"]) == "sinlike") {
            if (!(basis.contains("parameter") && basis["parameter"].is_number())) {
                MustContain("parameter", "number", "basis");
                return false;
            }
        } else if (ToLower(basis["node_sequence"]) == "parabolic") {
            if (!(basis.contains("parameter") && basis["parameter"].is_number())) {
                MustContain("parameter", "number", "basis");
                return false;
            }
        } else {
            LOG_CRITICAL("node_sequence not supported!");
            return false;
        }
        // CHECK SEQUENCES
        if (!(basis.contains("num_nodes") && basis["num_nodes"].is_number())) {
            MustContain("num_nodes", "number");
            return false;
        }
        if (!(basis.contains("x_min") && basis["x_min"].is_number())) {
            MustContain("x_min", "number");
            return false;
        }
        if (!(basis.contains("x_max") && basis["x_max"].is_number())) {
            MustContain("x_max", "number");
            return false;
        }
        if (!(basis.contains("lmax") && basis["lmax"].is_number())) {
            MustContain("lmax", "number");
            return false;
        }
        if (!(basis.contains("mmax") && basis["mmax"].is_number())) {
            MustContain("mmax", "number");
            return false;
        }
        if (basis["mmax"] > basis["lmax"]) {
            Log::critical("mmax must be less than or equal to lmax.");
            return false;
        }
    }
    return true;
}