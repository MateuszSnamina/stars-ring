#include<cmath>

#include<maths_in_physic/oscilator_realm.hpp>

namespace {

    arma::vec create_n_diag(unsigned n_max_stars) {
        arma::vec n_diag(n_max_stars);
        for (unsigned i = 0; i < n_max_stars; ++i)
            n_diag(i) = i;
        return n_diag;
    }

    arma::vec create_cr_an_diag(unsigned n_max_stars) {
        arma::vec cr_an_diag(n_max_stars - 1);
        for (unsigned i = 0; i < n_max_stars - 1; ++i)
            cr_an_diag(i) = sqrt(i + 1.0);
        return cr_an_diag;
    }

}

namespace maths_in_physic {

    OscilatorRealm::OscilatorRealm(unsigned requested_n_max_stars) :
    n_max_stars(requested_n_max_stars),
    n_diag(create_n_diag(n_max_stars)),
    cr_an_diag(create_cr_an_diag(n_max_stars)),
    n(arma::diagmat(n_diag)),
    cr(arma::diagmat(cr_an_diag, -1)),
    an(arma::diagmat(cr_an_diag, +1)) {
    }

    std::ostream& operator<<(std::ostream& s, const OscilatorRealm& oscilator_realm) {
        s << "OscilatorRealm::n:" << std::endl << oscilator_realm.n << std::endl;
        s << "OscilatorRealm::cr:" << std::endl << oscilator_realm.cr << std::endl;
        s << "OscilatorRealm::an:" << std::endl << oscilator_realm.an << std::endl;
        return s;
    }

}