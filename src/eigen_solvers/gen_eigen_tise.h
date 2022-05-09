#pragma once

#include "tise/tise.h"

class GeneralizeEigenvalueTISE : public TISE {
    MathLib& _MathLib;

    std::vector<std::vector<Vector>> _vectors;
    std::vector<std::vector<complex>> _values;
public:
    GeneralizeEigenvalueTISE(MathLib& mathlib);
    
    void Solve();
    void Store() const;

    void Finish();
};


