#include "input_validation/validate.h"

bool ValidateMathLibrary(const nlohmann::json& input) {
    if (!(input.contains("math_library") && input["math_library"].is_string())) {
        MustContain("math_library", "string");
        return false;
    }
    auto mathlibstr = ToLower(input["math_library"]);
    if (mathlibstr != "petsc" && mathlibstr != "threadpool") {
        return false;
    }
    return true;
}