#pragma once

#include <vector>
#include <memory>
#include "maths/maths.h"
#include "utility/vec3.h"

struct Sin2Pulse;

struct Pulse {
    typedef std::shared_ptr<Pulse> Ptr_t;

    const double num_cycles;
    const double E0, frequency, intensity;
    const double period, duration, delay, cep;
    const double ellipticity;
    const Vec3 polarization_vector;
    const Vec3 poynting_vector;
    const Vec3 minor_polarization_vector;
    const double cycles_plateau, cycles_up, cycles_down;

    Pulse();
    Pulse(double delay_cycles, double cep, double intensity, double frequency, double num_cycles, double cycles_up, double cycles_down, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);

    enum Envelope {
        Gaussian,
        Sin2,
        Box,
        Trap,
        
        NumTypes
    };

    virtual Vec3 operator() (double t) const = 0;
    virtual Vec3 A(double t) const = 0;
    virtual Vec3 E(double t) const = 0;

    static Ptr_t Create(Envelope env, double delay_cycles, double cep, double intensity, double frequency, double num_cycles, double cycles_up, double cycles_down, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);
};

struct Sin2Pulse : public Pulse {
    // probably add phase, time shift, etc.

    Sin2Pulse(double delay_cycles, double cep, double intensity, double frequency, 
             double num_cycles, double cycles_up, double cycles_down,
             double ellipticity,
             const Vec3& polarization,
             const Vec3& poynting_vector) :
        Pulse(delay_cycles, cep, intensity, frequency, num_cycles, cycles_up, cycles_down, ellipticity, polarization, poynting_vector) {}

    Vec3 operator() (double t) const;
    Vec3 A(double t) const;
    Vec3 E(double t) const;
};

struct TrapezoidalPulse : public Pulse {
    // probably add phase, time shift, etc.

    TrapezoidalPulse(
             double delay_cycles, double cep, double intensity, 
             double frequency, double num_cycles, double cycles_up, double cycles_down, 
             double ellipticity,
             const Vec3& polarization,
             const Vec3& poynting_vector) :
        Pulse(delay_cycles, cep, intensity, frequency, num_cycles, cycles_up, cycles_down, ellipticity, polarization, poynting_vector) {}

    Vec3 operator() (double t) const;
    Vec3 A(double t) const;
    Vec3 E(double t) const;
};