#include "input_validation/validate.h"

bool ValidatePropagator(const nlohmann::json& input) {
    if (!(input.contains("propagator") && input["propagator"].is_string())) {
        MustContain("propagator", "string");
        return false;
    }
    if (!(input.contains("time_step") && input["time_step"].is_number())) {
        MustContain("time_step", "number");
        return false;
    }
    return true;
}