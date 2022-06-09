#pragma once

#include <json.hpp>
#include <string>
#include "utility/logger.h"
#include "maths/maths.h"
#include "tdse/tdse.h"
#include "tise/tise.h"
#include "tdse/observable.h"
#include "tdse/potential.h"

inline void MustContain(const std::string& entry, const std::string& type = "", const std::string& location = "input file") {
    LOG_CRITICAL(location + " must contain " + type + " entry: " + entry);
}

inline std::string ToLower(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(),
        [](unsigned char c){ return std::tolower(c); });
    return data;
}

bool ValidatePotentials(const nlohmann::json& input);
bool ValidateEigenStateCalculation(const nlohmann::json& input);
bool ValidatePropagator(const nlohmann::json& input);
bool ValidateLasers(const nlohmann::json& input);
bool ValidateInitialState(const nlohmann::json& input);
bool ValidateMathLibrary(const nlohmann::json& input);
bool ValidateBasis(const nlohmann::json& input);
bool ValidateObservables(const nlohmann::json& input);

bool ValidateTISEInputFile(int argc, char **args, const std::string& filename, MathLib*& matlib, TISE::Ptr_t& tise);
bool ValidateTDSEInputFile(int argc, char **args, const std::string& filename, MathLib*& matlib, TDSE::Ptr_t& tdse);


Observable::Ptr_t BuildObservable(const std::string& key, const nlohmann::json& obs_item, TDSE::Ptr_t tdse);
Potential::Ptr_t BuildPotential(const nlohmann::json& potential_item);
