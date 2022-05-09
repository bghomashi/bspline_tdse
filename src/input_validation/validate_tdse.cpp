#include "input_validation/validate.h"

#include <iostream>
#include <fstream>

// math libraries
#include "math_libs/petsc/petsc_lib.h"

// propagators
#include "tdse_propagators/cranknicolson.h"

// observables
#include "observables/norm_obs.h"
#include "observables/dipole_acc_obs.h"

// potentials

bool ValidateTDSEInputFile(int argc, char **args, const std::string& filename, MathLib*& matlib, TDSE::Ptr_t& tdse) {
    std::ifstream i(filename);
    nlohmann::json input;
    i >> input;

    // log filename is optional. If not found output to stdout
    if (input.contains("log_filename") && input["log_filename"].is_string())
        Log::set_logger_file(input["log_filename"]);
        
    // first we immediately check the math library
    // - if this is library has a special 'logger' we start it right now
    if (!ValidateMathLibrary(input))
        return false;

    if (ToLower(input["math_library"]) == "petsc") {
        matlib = &Petsc::get();
        Log::set_logger(new PetscLogger());
    } else if (ToLower(input["math_library"]) == "thread_pool") {
        std::cout << "thread_pool is not yet supported" << std::endl;
        return false;
    }
    matlib->Startup(argc, args);


    Log::info("Validating TDSE input file.");

    Log::info("...");
    
    if (!ValidatePropagator(input))
        return false;
    Log::info("...");

    if (!ValidateBasis(input))
        return false;
    Log::info("...");

    if (!ValidateLasers(input))
        return false;
    Log::info("...");

    if (!ValidateInitialState(input))
        return false;
    Log::info("...");

    if (!ValidateObservables(input))
        return false;
    Log::info("...");
    if (!ValidatePotentials(input))
        return false;
    Log::info("...");
        

    // only need the "filename" entry from eigenstate calculation
    // assumes it is a valid filename since eigenstate calculation went through
    if (!(input.contains("eigen_state") && input["eigen_state"].is_object())) {
        MustContain("eigen_state", "object");
        return false;
    }
    if (!(input["eigen_state"].contains("filename") && input["eigen_state"]["filename"].is_string())) {
        MustContain("filename", "string");
        return false;
    }
    if (!(input.contains("time_step") && input["time_step"].is_number())) {
        MustContain("time_step", "number");
        return false;
    }
    if (!(input.contains("checkpoint") && input["checkpoint"].is_number())) {
        MustContain("checkpoint", "number");
        return false;
    }
    Log::info("...");

    Log::info("Input file validated. Initializing TDSE.");
    // we made it! begin loading data from input file
    // ORDER MATTERS - basis needs to know about lasers and initial states to 
    // take advatage of symmetry
    

    if (ToLower(input["propagator"]) == "crank_nicolson")
        tdse = TDSE::Ptr_t(new CrankNicolsonTDSE(*matlib));
    tdse->SetTimestep(input["time_step"]);
    tdse->SetCheckpoints(input["checkpoint"]);
    // set up lasers
    Log::info("Setting up lasers.");

    auto& lasers = input["lasers"];
    for (auto& pulse : lasers) {
        double frequency;
        double num_cycles = pulse["num_cycles"];
        double cycles_delay = pulse["cycles_delay"];
        double intensity = pulse["intensity"];
        double cep = pulse["cep"];
        std::vector<double> pol_vector(3);
        double norm = 0;
        for (int i = 0; i < 3; i++) {
            pol_vector[i] = pulse["polarization_vector"][i];
            norm += pol_vector[i]*pol_vector[i];
        }
        // not sure if we need to bother with normalization
        for (int i = 0; i < 3; i++) 
            pol_vector[i] /= sqrt(norm);

        // specify wavelength (nm) or energy (au)
        if (pulse.contains("wavelength"))
            frequency = LnmToEnergy/pulse["wavelength"].get<double>();   // not sure why "get..." is needed here
        if (pulse.contains("energy"))
            frequency = pulse["energy"];


        if (pulse["envelope"] == "sin2")
            tdse->AddPulse(Pulse::Create(Pulse::Sin2, cycles_delay, intensity, frequency, num_cycles, pol_vector));
    }
    
    // setup initial state
    Log::info("Setting up initial state.");


    tdse->SetInitialStateFile(input["eigen_state"]["filename"]);
    tdse->SetEigenStateNmax(input["eigen_state"]["nmax"]);
    auto& initial_state = input["initial_state"];
    for (auto& state : initial_state) {
        int n = state["n"];
        int l = state["l"];
        int m = state["m"];
        double phase = 0.;
        double amplitude = 1.;

        if (state.contains("phase")) phase = state["phase"];
        if (state.contains("amplitude")) amplitude = state["amplitude"];
        
        tdse->AddInitialState(n, l, m, phase, amplitude);
    }


    // setup basis
    Log::info("Setting up basis.");

    auto& basis = input["basis"];
    int order = basis["order"];
    int num_nodes = basis["num_nodes"];
    double x_min = basis["x_min"];
    double x_max = basis["x_max"];
    int lmax = basis["lmax"];
    int mmax = basis["mmax"];
    double ecs_r0 = 0.9;
    double ecs_theta = Pi/4.;
    Basis::BSpline::Sequence seq = Basis::BSpline::Linear;


    if (basis.contains("ecs_r0")) ecs_r0 = basis["ecs_r0"];           // optional - ecs parameters
    if (basis.contains("ecs_theta")) ecs_theta = basis["ecs_theta"];
    
    if (ToLower(basis["node_sequence"]) == "linear")
        seq = Basis::BSpline::Linear;
    else if (ToLower(basis["node_sequence"]) == "exponential")
        seq = Basis::BSpline::Exponential;
    else if (ToLower(basis["node_sequence"]) == "sinlike")
        seq = Basis::BSpline::Sinlike;
    else if (ToLower(basis["node_sequence"]) == "parabolic")
        seq = Basis::BSpline::ParabolicLinear;

    tdse->SetupBasis(x_min, x_max, 
                    order, num_nodes, 
                    lmax, mmax,
                    ecs_r0, ecs_theta,
                    seq);
    
    // setup observables
    Log::info("Building observables.");

    auto& observables_json = input["observables"];
    for (auto& obs_pair : observables_json.items()) {
        Observable::Ptr_t obs_ptr;
        if ((obs_ptr = BuildObservable(obs_pair.key(), obs_pair.value(), tdse)) == nullptr)
            return false;                                       // should never happen because we validated
        tdse->AddObservable(obs_ptr);
    }

    // setup potentials
    Log::info("Building potentials.");
    if (input.contains("potentials")) {
        auto& potentials = input["potentials"];
        for (auto& term : potentials) {
            Potential::Ptr_t pot_ptr;
            if ((pot_ptr = BuildPotential(term)) == nullptr)
                return false;                                       // should never happen because we validated
            tdse->AddPotential(pot_ptr);
        }
    }

    Log::info("Success.\n\n");
    return true;
}