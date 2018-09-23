#ifndef OSCILATOR_REALM_HPP
#define OSCILATOR_REALM_HPP

#include <armadillo>
#include <iostream>

namespace maths_in_physic {

struct OscilatorRealm {
  explicit OscilatorRealm(unsigned requested_n_max_stars);
  const unsigned n_max_stars;
  const arma::vec n_diag;
  const arma::vec cr_an_diag;
  const arma::mat n;
  const arma::mat cr;
  const arma::mat an;
};

std::ostream& operator<<(std::ostream&, const OscilatorRealm&);

}  // namespace maths_in_physic

#endif