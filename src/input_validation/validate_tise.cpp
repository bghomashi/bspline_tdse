
#include "input_validation/validate.h"
#include <iostream>
#include <fstream>
#include "utility/logger.h"
#include "utility/profiler.h"
#include "math_libs/petsc/petsc_lib.h"
#include "eigen_solvers/gen_eigen_tise.h"

bool ValidateTISEInputFile(int argc, char **args, const std::string& filename, MathLib*& matlib, TISE::Ptr_t& tise) {
    std::ifstream i(filename);
    nlohmann::json input;
    i >> input;

    // log filename is optional. If not found output to stdout
    if (input.contains("log_filename") && input["log_filename"].is_string())
        Log::set_logger_file(input["log_filename"]);

    auto& eigen_state = input["eigen_state"];
    if (!(eigen_state.contains("solver") && eigen_state["solver"].is_string())) {
        MustContain("solver", "string");
        return false;
    } else {
        if (ToLower(eigen_state["solver"]) == "slepc") {
            matlib = &Petsc::get();
            Log::set_logger(new PetscLogger());
            Profile::SetProfiler(new PetscProfiler());
        } else { //if (ToLower(eigen_state["solver"]) == "??") {
            Log::critical("only SLEPC is supported");
            return false;
        }
        matlib->Startup(argc, args);
    }


    Log::info("Validating TISE input file.");

    if (!ValidateEigenStateCalculation(input))
        return false;
    Log::info("...");
    if (!ValidateBasis(input))
        return false;
    Log::info("...");


    Log::info("Input file validated. Initializing TISE.");

    tise = TISE::Ptr_t(new GeneralizeEigenvalueTISE(*matlib));
    if (eigen_state.contains("tol"))
        tise->SetTolerance(eigen_state["tol"]);            // optional
    tise->SetNMax(eigen_state["nmax"]);
    tise->SetFilename(eigen_state["filename"]);
    // setup basis
    Log::info("Setting up basis.");

    auto& basis = input["basis"];
    int order = basis["order"];
    int num_nodes = basis["num_nodes"];
    double x_min = basis["x_min"];
    double x_max = basis["x_max"];
    int lmax = basis["lmax"];
    int mmax = basis["mmax"];
    double ecs_r0 = 1.0;                // default to NO ecs
    double ecs_theta = 0.0;             // default to NO ecs
    Basis::BSpline::Sequence seq = Basis::BSpline::Linear;  // default to linear

    // optional - ecs parameters
    if (basis.contains("ecs_r0")) ecs_r0 = basis["ecs_r0"];
    if (basis.contains("ecs_theta")) ecs_theta = basis["ecs_theta"];
    
    if (ToLower(basis["node_sequence"]) == "linear")
        seq = Basis::BSpline::Linear;
    else if (ToLower(basis["node_sequence"]) == "exponential")
        seq = Basis::BSpline::Exponential;
    else if (ToLower(basis["node_sequence"]) == "sinlike")
        seq = Basis::BSpline::Sinlike;
    else if (ToLower(basis["node_sequence"]) == "parabolic")
        seq = Basis::BSpline::ParabolicLinear;

    tise->SetupBasis(x_min, x_max, 
                    order, num_nodes, 
                    ecs_r0, ecs_theta,
                    seq);
    
    // setup potentials
    Log::info("Building potentials.");
    if (input.contains("potentials")) {
        auto& potentials = input["potentials"];
        for (auto& term : potentials) {
            Potential::Ptr_t pot_ptr;
            if ((pot_ptr = BuildPotential(term)) == nullptr)
                return false;                                       // should never happen because we validated
            tise->AddPotential(pot_ptr);
        }
    }

    Log::info("Success.\n\n");
    return true;
}