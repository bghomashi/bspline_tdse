#pragma once

#include <vector>
#include <memory>
#include "maths/maths.h"
#include "utility/vec3.h"

struct Sin2Pulse;

struct Pulse {
    typedef std::shared_ptr<Pulse> Ptr_t;

    const int numCycles;
    const double E0, frequency, intensity;
    const double period, duration, delay, cep;
    const double ellipticity;
    const Vec3 polarization_vector;
    const Vec3 poynting_vector;
    const Vec3 minor_polarization_vector;

    Pulse();
    Pulse(double delay_cycles, double cep, double intensity, double frequency, int numCycles, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);

    enum Envelope {
        Gaussian,
        Sin2,
        Box,
        
        NumTypes
    };

    virtual Vec3 operator() (double t) const = 0;
    virtual Vec3 A(double t) const = 0;
    virtual Vec3 E(double t) const = 0;

    static Ptr_t Create(Envelope env, double delay_cycles, double cep, double intensity, double frequency, int numCycles, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);
};

struct Sin2Pulse : public Pulse {
    // probably add phase, time shift, etc.

    Sin2Pulse(double delay_cycles, double cep, double intensity, 
             double frequency, int numCycles, 
             double ellipticity,
             const Vec3& polarization,
             const Vec3& poynting_vector) :
        Pulse(delay_cycles, cep, intensity, frequency, numCycles, ellipticity, polarization, poynting_vector) {}

    Vec3 operator() (double t) const;
    Vec3 A(double t) const;
    Vec3 E(double t) const;
};