
#include "tdse/laser.h"


Pulse::Pulse() : period(0), delay(0), numCycles(0), E0(0), duration(0), frequency(0), intensity(0), ellipticity(0), 
polarization_vector(), poynting_vector(), minor_polarization_vector()
{}
Pulse::Pulse(double delay_cycles, double intensity, 
             double frequency, int numCycles, 
             double ellipticity,
             const Vec3& polarization,
             const Vec3& poynting_vector) : 
    numCycles(numCycles), frequency(frequency), 
    intensity(intensity),
    E0(std::sqrt(intensity / 3.51e16)), 
    period(2.*Pi/frequency),
    delay(period*delay_cycles),
    duration(numCycles*period),
    ellipticity(ellipticity),
    polarization_vector(normal(polarization)/std::sqrt(1.+ellipticity*ellipticity)),
    poynting_vector(normal(poynting_vector)),
    minor_polarization_vector(normal(cross(polarization, poynting_vector))*(ellipticity/std::sqrt(1.+ellipticity*ellipticity)))
    { }
    
Pulse::Ptr_t Pulse::Create(Envelope env, 
                           double delay_cycles, double intensity, 
                           double frequency, int numCycles, 
                           double ellipticity,
                           const Vec3& polarization,
                           const Vec3& poynting_vector) {
    switch (env)
    {
    case Gaussian:
        // not yet supported
        break;
    case Sin2:
        return Pulse::Ptr_t(new Sin2Pulse(delay_cycles, intensity, frequency, numCycles, ellipticity, polarization, poynting_vector));
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

Vec3 Sin2Pulse::operator() (double t) {
    if (t < delay) return Vec3{0};

    // also need to add CEP
    double T = t-delay;
    double Env = E0*sin(Pi*T/duration)*sin(Pi*T/duration)/frequency;
    // Vec3 p = cos(frequency*T)*polarization_vector - sin(frequency*T)*minor_polarization_vector;
    // Vec3 p = sin(frequency*T)*polarization_vector - cos(frequency*T)*minor_polarization_vector;      // try this next
    Vec3 p = sin(frequency*T)*polarization_vector;
    return Env*p;
}