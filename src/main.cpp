#include "tdse_propagators/cranknicolson.h"
#include "eigen_solvers/gen_eigen_tise.h"
#include "math_libs/petsc/petsc_lib.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <json.hpp>
#include "input_validation/validate.h"
#include "core/utility/logger.h"


#include "core/utility/profiler.h"

int main(int argc, char **args) {
    MathLib* matLib;
    bool do_tise = false;
    TDSE::Ptr_t tdse;
    TISE::Ptr_t tise;

    // solve TISE?
    for (int i = 0; i < argc; i++) {
        if (ToLower(std::string(args[i])) == "tise" ||
            ToLower(std::string(args[i])) == "-tise" ||
            ToLower(std::string(args[i])) == "--tise")
            do_tise = true;
    }

    if (do_tise) {
        Profile::Push("Total TISE time");
        if (!ValidateTISEInputFile(argc, args, "input.json", matLib, tise))
            return -1;

        LOG_INFO("Starting up TISE");
        tise->Solve();
        tise->Finish();
        Profile::Pop("Total TISE time");
    } else {
        Profile::Push("Total TDSE time");
        if (!ValidateTDSEInputFile(argc, args, "input.json", matLib, tdse))
            return -1;
        
        LOG_INFO("Starting up TDSE");
        tdse->Initialize();
        tdse->Propagate();
        Profile::Pop("Total TDSE time");
    }
    LOG_INFO("Shutting down.\n------------------------------------------------\n\n");
    Profile::PrintTo("profile.txt");
    
    matLib->Shutdown();

    return 0;
}