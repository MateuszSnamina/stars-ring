#include <armadillo>
#include <array>
#include <cmath>

#include <stars_ring_analytical/analytical_formulas_box_af_spins.hpp>
#include <stars_ring_analytical/analytical_formulas_box_jabcdz.hpp>

namespace {

struct ArgAndExpansion {
  double x, f, fp, fpp;  // values x, f(x), f'(x), f''(x)
};

class AlphaCos2PlusBetaCos {
 public:
  AlphaCos2PlusBetaCos(double alpha, double beta);
  double f(double x) const;
  double fp(double x) const;
  double fpp(double x) const;
  ArgAndExpansion arg_and_expansion(double x) const;
  ArgAndExpansion get_minimum() const;

 private:
  const double _alpha, _beta;
  const std::array<ArgAndExpansion, 3>
      _data;  // three values for which f'(x)=0;
};

AlphaCos2PlusBetaCos::AlphaCos2PlusBetaCos(double alpha, double beta)
    : _alpha(alpha),
      _beta(beta),
      _data{arg_and_expansion(0.0),
            arg_and_expansion(std::acos(-_beta / 2.0 / _alpha)),
            arg_and_expansion(arma::datum::pi)} {
  std::cout << "[DEBUG  ] minimizer alpha, beta   : " << _alpha << " " << _beta
            << "." << std::endl;
  std::cout << "[DEBUG  ] minimizer x, f, f', f'' : " << _data[0].x << " "
            << _data[0].f << " " << _data[0].fp << " " << _data[0].fpp << "."
            << std::endl;
  std::cout << "[DEBUG  ] minimizer x, f, f', f'' : " << _data[1].x << " "
            << _data[1].f << " " << _data[1].fp << " " << _data[1].fpp << "."
            << std::endl;
  std::cout << "[DEBUG  ] minimizer x, f, f', f'' : " << _data[2].x << " "
            << _data[2].f << " " << _data[2].fp << " " << _data[2].fpp << "."
            << std::endl;
}

double AlphaCos2PlusBetaCos::f(double x) const {
  return _alpha * std::pow(std::cos(x), 2) + _beta * std::cos(x);
}

double AlphaCos2PlusBetaCos::fp(double x) const {
  return -(2 * _alpha * std::cos(x) + _beta) * std::sin(x);
}

double AlphaCos2PlusBetaCos::fpp(double x) const {
  return 2 * _alpha * std::pow(std::sin(x), 2) -
         (2 * _alpha * std::cos(x) + _beta) * std::cos(x);
}

ArgAndExpansion AlphaCos2PlusBetaCos::arg_and_expansion(double x) const {
  return {x, f(x), fp(x), fpp(x)};
}

ArgAndExpansion AlphaCos2PlusBetaCos::get_minimum() const {
  if (_alpha == 0.0 && _beta == 0.0) {
    std::cerr << "[WARNING] cannnot mimimize f(x) = 0 * cos^2(x) + 0 * cos(x) "
                 "function."
              << std::endl;
  }
  ArgAndExpansion result;
  result = _data[0];
  for (unsigned i = 1; i < _data.size(); i++) {
    if (result.f > _data[i].f) result = _data[i];
  }
  assert(!std::isnan(result.x));
  assert(std::abs(result.fp) < 1e-7);
  assert(result.fpp > -1e-7);
  return result;
}
/*
    double calculate_theta_opt_nell(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            double multiplicity,
            double A, double B, double D,
            double J, double Ez) {
        const stars_ring_analytical::AnalyticalFormulasBoxAfSpins
   analytical_formulas_box_af_spins(physical_system, multiplicity); const double
   alpha = J * A * D * physical_system->n_sites() + J * B * D *
   analytical_formulas_box_af_spins.ground_state_classical_energy(); const
   double beta = physical_system->n_sites() * (-Ez) * (-1.0 / 2.0);
        AlphaCos2PlusBetaCos fun(alpha, beta);
        return fun.get_minimum().x;
    }

    double calculate_theta_opt_corrected_nell(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            double multiplicity,
            double A, double B, double D,
            double J, double Ez) {
        const stars_ring_analytical::AnalyticalFormulasBoxAfSpins
   analytical_formulas_box_af_spins(physical_system, multiplicity); const double
   alpha = J * A * D * physical_system->n_sites() + J * B * D *
   analytical_formulas_box_af_spins.ground_state_energy(); const double beta =
   physical_system->n_sites() * (-Ez) * (-1.0 / 2.0); AlphaCos2PlusBetaCos
   fun(alpha, beta); return fun.get_minimum().x;
    }
 */
}  // namespace

namespace stars_ring_analytical {

AnalyticalFormulasBoxJABCDZ::AnalyticalFormulasBoxJABCDZ(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    unsigned multiplicity, double A, double B, double C, double D, double J,
    double Ez, double theta)
    : AnalyticalFormulasBox(physical_system),
      _analytical_formulas_box_af_spins(physical_system, multiplicity),
      _multiplicity(multiplicity),
      _A(A),
      _B(B),
      _C(C),
      _D(D),
      _J(J),
      _Ez(Ez),
      _theta(theta) {
  std::cout << "[DEBUG  ] [AnalyticalFormulasBoxSu4] A, B, C, D : " << _A
            << ", " << _B << ", " << _C << ", " << _D << std::endl;
  std::cout << "[DEBUG  ] [AnalyticalFormulasBoxSu4] J, Ez      : " << _J
            << ", " << _Ez << std::endl;
}

// ponizsza funkcja powinna byc funkcją konstruktorowa, nie konstruktorem:
/*
AnalyticalFormulasBoxJABCDZ::AnalyticalFormulasBoxJABCDZ(
        std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
        unsigned multiplicity,
        double A, double B, double C, double D,
        double J, double Ez) :
AnalyticalFormulasBoxJABCDZ(
physical_system,
multiplicity,
A, B, C, D,
J, Ez,
calculate_theta_opt_total(physical_system, multiplicity, A, B, D, J, Ez)) {
}
 */

double AnalyticalFormulasBoxJABCDZ::ground_state_classical_energy() const {
  return _J * _A * mean_orbital_operator() * physical_system()->n_sites() +
         J_spin() *
             _analytical_formulas_box_af_spins.ground_state_classical_energy() -
         _Ez * (-std::cos(_theta) / 2.0) * physical_system()->n_sites();
}

double AnalyticalFormulasBoxJABCDZ::ground_state_correlation_energy() const {
  return J_spin() *
         _analytical_formulas_box_af_spins.ground_state_correlation_energy();
}

double AnalyticalFormulasBoxJABCDZ::exc_state_relative_energy(
    unsigned nk) const {
  return J_spin() *
         _analytical_formulas_box_af_spins.exc_state_relative_energy(nk);
}

double AnalyticalFormulasBoxJABCDZ::theta() const { return _theta; }

double AnalyticalFormulasBoxJABCDZ::mean_orbital_operator() const {
  return _C + _D * std::pow(std::cos(_theta), 2) / 4.0;
}

double AnalyticalFormulasBoxJABCDZ::J_spin() const {
  return _J * _B * mean_orbital_operator();
}

double AnalyticalFormulasBoxJABCDZ::theta_opt_nell() const {
  const stars_ring_analytical::AnalyticalFormulasBoxAfSpins
      analytical_formulas_box_af_spins(physical_system(), _multiplicity);
  const double alpha =
      _J * _A * _D * physical_system()->n_sites() +
      _J * _B * _D *
          analytical_formulas_box_af_spins.ground_state_classical_energy();
  const double beta = physical_system()->n_sites() * (-_Ez) * (-1.0 / 2.0);
  AlphaCos2PlusBetaCos fun(alpha, beta);
  return fun.get_minimum().x;
}

double AnalyticalFormulasBoxJABCDZ::theta_opt_corrected_nell() const {
  const stars_ring_analytical::AnalyticalFormulasBoxAfSpins
      analytical_formulas_box_af_spins(physical_system(), _multiplicity);
  const double alpha =
      _J * _A * _D * physical_system()->n_sites() +
      _J * _B * _D * analytical_formulas_box_af_spins.ground_state_energy();
  const double beta = physical_system()->n_sites() * (-_Ez) * (-1.0 / 2.0);
  AlphaCos2PlusBetaCos fun(alpha, beta);
  return fun.get_minimum().x;
}

}  // namespace stars_ring_analytical