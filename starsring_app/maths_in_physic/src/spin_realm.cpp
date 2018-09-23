#include <cmath>

#include <maths_in_physic/spin_realm.hpp>

namespace {

arma::vec create_s_z_diag(unsigned multiplicity) {
  const unsigned doubled_l = multiplicity - 1;
  arma::vec S_z_diag(multiplicity);
  for (unsigned n_stars = 0; n_stars < multiplicity; n_stars++) {
    const int doubled_m = int(doubled_l) - 2 * n_stars;
    S_z_diag(n_stars) = doubled_m / 2.0;
  }
  return S_z_diag;
}

arma::vec create_s_pm_diag(unsigned multiplicity) {
  const unsigned doubled_l = multiplicity - 1;
  arma::vec S_pm_diag(multiplicity - 1);
  for (unsigned n_stars = 0; n_stars < multiplicity - 1; n_stars++) {
    const int doubled_m = int(doubled_l) - 2 * n_stars;
    S_pm_diag(n_stars) =
        std::sqrt(doubled_l * (doubled_l + 2) - (doubled_m - 2) * doubled_m) /
        2.0;
  }
  return S_pm_diag;
}

}  // namespace

namespace maths_in_physic {

SpinRealm::SpinRealm(unsigned requested_multiplicity)
    : multiplicity(requested_multiplicity),
      doubledS(multiplicity - 1),
      S(doubledS / 2.0),
      S_z_diag(create_s_z_diag(multiplicity)),
      S_pm_diag(create_s_pm_diag(multiplicity)),
      S_z(arma::diagmat(S_z_diag)),
      S_p(arma::diagmat(S_pm_diag, +1)),
      S_m(arma::diagmat(S_pm_diag, -1)) {}

std::ostream& operator<<(std::ostream& s, const SpinRealm& spin_realm) {
  s << "S_z:" << std::endl << spin_realm.S_z << std::endl;
  s << "S^+:" << std::endl << spin_realm.S_p << std::endl;
  s << "S^-:" << std::endl << spin_realm.S_m << std::endl;
  return s;
}

}  // namespace maths_in_physic