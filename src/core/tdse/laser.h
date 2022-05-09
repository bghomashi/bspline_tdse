#pragma once

#include <vector>
#include <memory>
#include "maths/maths.h"

struct Sin2Pulse;

struct Pulse {
    typedef std::shared_ptr<Pulse> Ptr_t;

    const int numCycles;
    const double E0, frequency, intensity;
    const double period, duration, delay;
    const std::vector<double> polarization_vector;

    virtual double operator() (double t) = 0;

    Pulse();
    Pulse(double delay_cycles, double intensity, double frequency, int numCycles, const std::vector<double>& polarization);

    enum Envelope {
        Gaussian,
        Sin2,
        Box,
        
        NumTypes
    };

    static Ptr_t Create(Envelope env, double delay_cycles, double intensity, double frequency, int numCycles, const std::vector<double>& polarization);
};

struct Sin2Pulse : public Pulse {
    // probably add phase, time shift, etc.

    Sin2Pulse(double delay_cycles, double intensity, double frequency, int numCycles, const std::vector<double>& polarization) : 
        Pulse(delay_cycles, intensity, frequency, numCycles, polarization) {}

    double operator() (double t);

};