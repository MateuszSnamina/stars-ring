#ifndef HAMILTONIAN_CREATOR_JABCDZ_HPP
#define HAMILTONIAN_CREATOR_JABCDZ_HPP

#include <armadillo>
#include <memory>

#include <maths_in_physic/spin_realm.hpp>
#include <stars_ring_basis/basis.hpp>
#include <stars_ring_numerical/phis_establisher.hpp>
#include <stars_ring_numerical/standard_calculator.hpp>

/*
 * Rozwazamy tu Hamiltonian:
 * H = J * (A + B * \vec \hat S_i \cdot \vec \hat S_j) *
 *         (C + D * \hat tau^c_i * \hat tau^c_j)
 *     - Ez * tau^c_i
 *
 * promujący AF/FO gdy Ez=0 (chyba...?)
 *
 * TCF:
 * Ez>0 promuje xi_c (theta = 180).
 * Ez<0 promuje zeta_c (theta = 0).
 *
 * Operatory orbitalne:
 * tau_c = ( |xi_c><xi_c| - |zeta_c><zeta_c| ) / 2.0   // Konwencja ,,po
 * staremu'' [przed zamiana tau := - tau] np. <zeta|tau_c|zeta> = <0|tau_c|0> =
 * -1 np. <xi|tau_c|xi>     = <180|tau_c|180> = +1
 *
 * Uwaga: dla C = 0.25, D = 1.0
 *        operator orbitalny to nie jest ani P^{xi,xi}, ani P^{zeta,xi}, ani
 * P^{zeta,zeta}, operaotr orbitalny P^{zeta,xi} to jest (1/4 - \hat tau^c_i *
 * \hat tau^c_j).
 *
 * Uawga: Dla S=1/2:
 *        P(tryplet) <=> A = +3/4, B = +1.0
 *        P(singlet) <=> A = -1/4, B = -1.0
 *
 * Uwaga: dla A = 0.25, B = 1.0
 *        operator spinowy to nie jest żaden z tych zwykłych.
 *        Przez ,,zwykle'' rozumiem operatory,
 *        Która dla S = 1/2 wyglądają następująco:
 *        P(tryplet) = (S*(S+1) + \vec \hat S \cdot \vec \hat S)
 *        P(singlet) = (S*S     - \vec \hat S \cdot \vec \hat S).
 *
 *
 * Klasa buduje hamiltonian po transformacji:
 * H = J * (A - B * \hat K^z_i * \hat K^z_j + 1/2 * B * \hat K^+_i * \hat K^+_j
 * + B * 1/2 * \hat K^-_i * \hat K^-_j) * (C + D * tau^c_i * tau^c_j)
 *     - Ez * tau^c_i
 * gdzie:
 * K_i = + S_i dla podsieci A, K_i = - S_i dla podsieci B.
 *
 * Dla modelu SU(4):
 * A = S^2
 * B = +1
 * C = 1/4
 * D = +1
 *
 */

namespace stars_ring_numerical {

class JABCDZHamiltonianCreator : public HamiltonianCreator {
 public:
  JABCDZHamiltonianCreator(
      maths_in_physic::SpinRealm spin_realm, double A, double B, double C,
      double D, double J, double Ez,
      std::shared_ptr<const PhisEstablisher> phis_establisher);
  virtual void creat_hamiltonian(
      arma::mat& H,
      const stars_ring_basis::LocalizedBasis& basis) const override;
  virtual void creat_hamiltonian(
      arma::sp_mat& H,
      const stars_ring_basis::LocalizedBasis& basis) const override;

 private:
  const maths_in_physic::SpinRealm _spin_realm;
  double _A, _B, _C, _D;
  double _J, _Ez;
  std::shared_ptr<const PhisEstablisher> _phis_establisher;
};

}  // namespace stars_ring_numerical

#endif
