#ifndef GLOBAL_SPIN_REALM_HPP
#define GLOBAL_SPIN_REALM_HPP

#include<iostream>
#include<armadillo>

namespace maths_in_physic {

    struct SpinRealm {
        explicit SpinRealm(unsigned requested_multiplicity);
        const unsigned multiplicity;
        const unsigned doubledS;
        const double S;
        const arma::vec S_z_diag;
        const arma::vec S_pm_diag;
        const arma::mat S_z;
        const arma::mat S_p;
        const arma::mat S_m;
    };

    std::ostream& operator<<(std::ostream&, const SpinRealm&);

}

#endif