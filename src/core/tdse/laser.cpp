
#include "tdse/laser.h"
#include "core/utility/logger.h"
#include <sstream>

Pulse::Pulse() : period(0), delay(0), numCycles(0), E0(0), duration(0), frequency(0), intensity(0), ellipticity(0), cep(0),
polarization_vector(), poynting_vector(), minor_polarization_vector()
{}
Pulse::Pulse(double delay_cycles, double cep, double intensity, 
             double frequency, int numCycles, 
             double ellipticity,
             const Vec3& polarization,
             const Vec3& poynting_vector) : 
    numCycles(numCycles), cep(cep*2.*Pi), frequency(frequency), 
    intensity(intensity),
    E0(std::sqrt(intensity / 3.51e16)), 
    period(2.*Pi/frequency),
    delay(period*delay_cycles),
    duration(numCycles*period),
    ellipticity(ellipticity),
    polarization_vector(normal(polarization)/std::sqrt(1.+ellipticity*ellipticity)),
    poynting_vector(normal(poynting_vector)),
    minor_polarization_vector(normal(cross(polarization, poynting_vector))*(ellipticity/std::sqrt(1.+ellipticity*ellipticity)))
    { 
        std::stringstream ss;

        ss << "Polarization = " << polarization_vector.x << ", " << polarization_vector.y << ", " << polarization_vector.z << std::endl;
        LOG_INFO(ss.str()); ss.str("");
        ss << "Poynting = " << poynting_vector.x << ", " << poynting_vector.y << ", " << poynting_vector.z << std::endl;
        LOG_INFO(ss.str()); ss.str("");
        ss << "MinorPolarization = " << minor_polarization_vector.x << ", " << minor_polarization_vector.y << ", " << minor_polarization_vector.z << std::endl;
        LOG_INFO(ss.str());
    }
    
Pulse::Ptr_t Pulse::Create(Envelope env, 
                           double delay_cycles, double cep,
                           double intensity, 
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
        return Pulse::Ptr_t(new Sin2Pulse(delay_cycles, cep, intensity, frequency, numCycles, ellipticity, polarization, poynting_vector));
        break;
    case Box:
        // not yet supported
        break;
    case Trap:
        return Pulse::Ptr_t(new TrapezoidalPulse(delay_cycles, cep, intensity, frequency, numCycles, ellipticity, polarization, poynting_vector));
        break;
    
    default:
        // LogError
        break;
    }
    return nullptr;
}

Vec3 Sin2Pulse::operator() (double t) const {
    return A(t);
}
Vec3 Sin2Pulse::A(double t) const {     // A(t)
    if (t < delay) return Vec3{0};

    double T = t-delay;
    double Env = -E0*sin(Pi*T/duration)*sin(Pi*T/duration)/frequency;
    Vec3 p = sin(frequency*T+ cep)*polarization_vector - cos(frequency*T + cep)*minor_polarization_vector;
    return Env*p;
}

Vec3 Sin2Pulse::E(double t) const {              // E=-dA/dt
    if (t < delay) return Vec3{0};

    double T = t-delay;
    // product rule
    double Env1 = E0*sin(Pi*T/duration)*sin(Pi*T/duration);                             // dont diff. env
    double Env2 = 2.*Pi*E0*sin(Pi*T/duration)*cos(Pi*T/duration)/frequency/duration;    // do diff. env

    Vec3 p1 = cos(frequency*T + cep)*polarization_vector + sin(frequency*T + cep)*minor_polarization_vector;       // do diff. carrier wave
    Vec3 p2 = sin(frequency*T + cep)*polarization_vector - cos(frequency*T + cep)*minor_polarization_vector;                           // dont

    return Env1*p1 + Env2*p2;

}





Vec3 TrapezoidalPulse::operator() (double t) const {
    return A(t);
}
Vec3 TrapezoidalPulse::A(double t) const {     // A(t)
    if (t < delay) return Vec3{0};

    auto envelope = [=](double t) {
        double t1 = cycles_up*2.*Pi/frequency;
        double t2 = cycles_plateau*2.*Pi/frequency + t1;
        double t3 = cycles_down*2.*Pi/frequency + t2;
        if (t <= 0)
            return 0.0;
        else if (t <= t1)
            return t/t1;
        else if (t <= t2)
            return 1.0;
        else if (t <= t3)
            return (t3-t) / (t3-t2);
        return 0.0;
    };

    double T = t-delay;
    double Env = -E0*envelope(T)/frequency;
    Vec3 p = sin(frequency*T+cep)*polarization_vector - cos(frequency*T+cep)*minor_polarization_vector;
    return Env*p;
}

Vec3 TrapezoidalPulse::E(double t) const {              // E=-dA/dt
    if (t < delay) return Vec3{0};

    auto envelope = [=](double t) {
        double t1 = cycles_up*2.*Pi/frequency;
        double t2 = cycles_plateau*2.*Pi/frequency + t1;
        double t3 = cycles_down*2.*Pi/frequency + t2;
        if (t <= 0)
            return 0.0;
        else if (t <= t1)
            return t/t1;
        else if (t <= t2)
            return 1.0;
        else if (t <= t3)
            return (t3-t) / (t3-t2);
        return 0.0;
    };
    auto denvelope = [=](double t) {
        double t1 = cycles_up*2.*Pi/frequency;
        double t2 = cycles_plateau*2.*Pi/frequency + t1;
        double t3 = cycles_down*2.*Pi/frequency + t2;
        if (t <= 0)
            return 0.0;
        else if (t <= t1)
            return 1./t1;
        else if (t <= t2)
            return 0.0;
        else if (t <= t3)
            return -1. / (t3-t2);
        return 0.0;
    };


    double T = t-delay;
    // product rule
    double Env1 = E0*envelope(T);                             // dont diff. env
    double Env2 = E0*denvelope(T)/frequency;    // do diff. env

    Vec3 p1 = cos(frequency*T+cep)*polarization_vector + sin(frequency*T+cep)*minor_polarization_vector;       // do diff. carrier wave
    Vec3 p2 = sin(frequency*T+cep)*polarization_vector - cos(frequency*T+cep)*minor_polarization_vector;                           // dont

    return Env1*p1 + Env2*p2;

}