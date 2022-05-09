#include "input_validation/validate.h"

bool ValidateInitialState(const nlohmann::json& input) {
    if (!(input.contains("initial_state") && input["initial_state"].is_array())) {
        MustContain("initial_state", "array");
        return false;
    }

    auto& initial_state = input["initial_state"];
    for (auto& state : initial_state) {
        if (!(state.contains("n") && state["n"].is_number())) {
            MustContain("n", "number");
            return false;
        }
        if (!(state.contains("l") && state["l"].is_number())) {
            MustContain("l", "number");
            return false;
        }
        if (!(state.contains("m") && state["m"].is_number())) {
            MustContain("m", "number");
            return false;
        }
        if (state.contains("phase") && !state["phase"].is_number()) {
            Log::critical("Optional entry \"phase\" must be a number.");
            return false;
        }
        if (state.contains("amplitude") && !state["amplitude"].is_number()) {
            Log::critical("Optional entry \"amplitude\" must be a number.");
            return false;
        }
    }
    return true;
}
