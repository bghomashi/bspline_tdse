#include "input_validation/validate.h"

bool ValidateEigenStateCalculation(const nlohmann::json& input) {
    if (!(input.contains("eigen_state") && input["eigen_state"].is_object())) {
        MustContain("eigen_state", "object");
        return false;
    }

    auto& eigen_state = input["eigen_state"];
    if (!(eigen_state.contains("solver") && eigen_state["solver"].is_string())) {
        MustContain("solver", "string");
        return false;
    }
    if (!(eigen_state.contains("nmax") && eigen_state["nmax"].is_number())) {
        MustContain("nmax", "number");
        return false;
    }
    if (!(eigen_state.contains("filename") && eigen_state["filename"].is_string())) {
        MustContain("filename", "string");
        return false;
    }
    if (eigen_state["nmax"] <= 0) {
        Log::critical("Entry nmax must be a number greater than 0.");
        return false;
    }
    if (eigen_state.contains("lmax") && !eigen_state["lmax"].is_number()) {
        Log::critical("Entry lmax must be a number greater than 0.");
        return false;
    }
    if (eigen_state.contains("lmin") && !eigen_state["lmin"].is_number()) {
        Log::critical("Entry lmin must be a number greater than 0.");
        return false;
    }
    {
        std::string filename = eigen_state["filename"];
        if ((filename.length() <= 3) || 
            (filename.find_last_of(".") == std::string::npos) ||
            (filename.substr(filename.find_last_of(".") + 1) != "h5")) {
            Log::critical("Filename entry must be at least one character *plus* \".h5\" extension");
            return false;
        }
    }

    return true;
}
