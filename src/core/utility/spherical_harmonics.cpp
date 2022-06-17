#include "utility/spherical_harmonics.h"

using namespace std::complex_literals;

complex YlmXYlm(int l1, int m1, int l2, int m2) {
    if (m1 == m2 + 1 && l1 == l2 + 1)
        return -0.5*sqrt( (l2+m2+1.)*(l2+m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 + 1 && l1 == l2 - 1)
        return 0.5*sqrt( (l2-m2)*(l2-m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    if (m1 == m2 - 1 && l1 == l2 + 1)
        return 0.5*sqrt( (l2-m2+1.)*(l2-m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 - 1 && l1 == l2 - 1)
        return -0.5*sqrt( (l2+m2)*(l2+m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    return 0.;
}
complex YlmYYlm(int l1, int m1, int l2, int m2) {                   // SIGN FLIP??
    if (m1 == m2 + 1 && l1 == l2 + 1)
        return 0.5i*sqrt( (l2+m2+1.)*(l2+m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 + 1 && l1 == l2 - 1)
        return -0.5i*sqrt( (l2-m2)*(l2-m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    if (m1 == m2 - 1 && l1 == l2 + 1)
        return 0.5i*sqrt( (l2-m2+1.)*(l2-m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 - 1 && l1 == l2 - 1)
        return -0.5i*sqrt( (l2+m2)*(l2+m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    return 0.;
}
complex YlmZYlm(int l1, int m1, int l2, int m2) {
    if (m1 != m2) return 0;
    if (l1 == l2 + 1) 
        return sqrt( (l2-m2+1.)*(l2+m2+1.)/(2.*l2+1.)/(2.*l2+3.) );
    
    if (l1 == l2 - 1) 
        return sqrt( (l2-m2)*(l2+m2)/(2.*l2+1)/(2.*l2-1) );
    

    return 0.;
}