
#include "tdse/laser.h"


Pulse::Pulse() : period(0), delay(0), numCycles(0), E0(0), duration(0), frequency(0), intensity(0) {}
Pulse::Pulse(double delay_cycles, double intensity, double frequency, int numCycles, const std::vector<double>& polarization) : 
    numCycles(numCycles), frequency(frequency), 
    intensity(intensity),
    E0(std::sqrt(intensity / 3.51e16)), 
    period(2.*Pi/frequency),
    delay(period*delay_cycles),
    duration(numCycles*period),
    polarization_vector(polarization) { }
    
Pulse::Ptr_t Pulse::Create(Envelope env, double delay_cycles, double intensity, double frequency, int numCycles, const std::vector<double>& polarization) {
    switch (env)
    {
    case Gaussian:
        // not yet supported
        break;
    case Sin2:
        return Pulse::Ptr_t(new Sin2Pulse(delay_cycles, intensity, frequency, numCycles, polarization));
        break;
    case Box:
        // not yet supported
        break;
    
    default:
        // LogError
        break;
    }
    return nullptr;
}

double Sin2Pulse::operator() (double t) {
    return E0*sin(Pi*t/duration)*sin(Pi*t/duration) * sin(frequency*t)/frequency;
}