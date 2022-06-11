#include "input_validation/validate.h"
#include "observables/norm_obs.h"
#include "observables/dipole_acc_obs.h"
#include "observables/wavefunction_obs.h"
#include "observables/population_obs.h"
#include "observables/pulse_obs.h"
#include "observables/potential_obs.h"
#include "observables/basis_obs.h"

Observable::Ptr_t BuildObservable(const std::string& key, const nlohmann::json& obs_item, TDSE::Ptr_t tdse) {
    // grab the compute period if provided
    int compute_period = 1;
    std::string filename;
    if (obs_item.contains("compute_period")) 
        compute_period = obs_item["compute_period"];
    if (obs_item.contains("filename")) 
        filename = obs_item["filename"];

    // make observable object
    if (key == "norm") {
        NormObservable* norm_obs = new NormObservable(*tdse);

        norm_obs->SetFilename(filename);
        norm_obs ->SetComputePeriod(compute_period);

        return Observable::Ptr_t(norm_obs);
    } else if (key == "dipole_acc") {
        DipoleAccObservable* dip_acc_obs = new DipoleAccObservable(*tdse);

        dip_acc_obs ->SetComputePeriod(compute_period);
        dip_acc_obs->SetFilename(filename);

        return Observable::Ptr_t(dip_acc_obs);
    } else if (key == "populations") {
        PopulationObservable* pop_obs = new PopulationObservable(*tdse);

        pop_obs->SetComputePeriod(compute_period);
        pop_obs->SetFilename(filename);

        return Observable::Ptr_t(pop_obs);
    } else if (key == "wavefunction") {
        WavefunctionObservable* wf_obs = new WavefunctionObservable(*tdse);

        wf_obs->SetComputePeriod(compute_period);
        wf_obs->SetFilename(filename);
        wf_obs->SetNumGrid(obs_item["grid_points"]);

        return Observable::Ptr_t(wf_obs);
    } else if (key == "pulse") {
        PulseObservable* pl_obs = new PulseObservable(*tdse);

        pl_obs->SetComputePeriod(compute_period);
        pl_obs->SetFilename(filename);

        return Observable::Ptr_t(pl_obs);
    } else if (key == "potential") {
        PotentialObservable* pot_obs = new PotentialObservable(*tdse);

        pot_obs->SetComputePeriod(compute_period);
        pot_obs->SetNumGrid(obs_item["grid_points"]);
        pot_obs->SetFilename(filename);

        return Observable::Ptr_t(pot_obs);
    } else if (key == "basis") {
        BasisObservable* basis_obs = new BasisObservable(*tdse);

        basis_obs->SetComputePeriod(compute_period);
        basis_obs->SetNumGrid(obs_item["grid_points"]);
        basis_obs->SetFilename(filename);

        return Observable::Ptr_t(basis_obs);
    }
    return nullptr;
}