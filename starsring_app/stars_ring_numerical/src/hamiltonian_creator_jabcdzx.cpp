#include <stdexcept>

#include <stars_ring_basis/raw_state.hpp>

#include <stars_ring_basis/raw_state_operations.hpp>
#include <stars_ring_numerical/hamiltonian_creator_jabcdzx.hpp>

namespace {

double phi1_id_psi2(double theta_1, double theta_2) {
  return +std::cos((theta_1 - theta_2) / 2.0);
}

double phi1_tauc_psi2(double theta_1, double theta_2) {
  // tau^c = + |xi><xi| - |zeta><zeta|
  return -std::cos((theta_1 + theta_2) / 2.0) / 2.0;
}

double phi1_tauanty_psi2(double theta_1, double theta_2) {
  // tau^{anty_c} = + |+><+| - |-><-|
  return +std::sin((theta_1 + theta_2) / 2.0) / 2.0;
}

double psi1_id_psi2(const RawPhis &bra_phis, const RawPhis &ket_phis) {
  assert(bra_phis.size() == ket_phis.size());
  double result = 1.0;
  for (unsigned idx = 0; idx < bra_phis.size(); ++idx)
    result *= phi1_id_psi2(bra_phis[idx], ket_phis[idx]);
  return result;
}

double psi1_tauci_taucj_psi2(const RawPhis &bra_phis, unsigned i, unsigned j,
                             const RawPhis &ket_phis) {
  assert(i != j);
  assert(i < bra_phis.size());
  assert(j < bra_phis.size());
  assert(bra_phis.size() == ket_phis.size());
  double result = 1.0;
  for (unsigned idx = 0; idx < bra_phis.size(); ++idx) {
    if (idx != i && idx != j)
      result *= phi1_id_psi2(bra_phis[idx], ket_phis[idx]);
    else
      result *= phi1_tauc_psi2(bra_phis[idx], ket_phis[idx]);
  }
  return result;
}

double psi1_orbitij_psi2(const RawPhis &bra_phis, const unsigned i,
                         const unsigned j, const RawPhis &ket_phis,
                         const double alpha, const double beta) {
  assert(i != j);
  assert(i < bra_phis.size());
  assert(j < bra_phis.size());
  assert(bra_phis.size() == ket_phis.size());
  double orbit_ss = psi1_id_psi2(bra_phis, ket_phis);
  double orbit_tt = psi1_tauci_taucj_psi2(bra_phis, i, j, ket_phis);
  return alpha * orbit_ss + beta * orbit_tt;
}

}  // namespace

namespace stars_ring_numerical {

JABCDZXHamiltonianCreator::JABCDZXHamiltonianCreator(
    maths_in_physic::SpinRealm spin_realm, double A, double B, double C,
    double D, double J, double Ez, double Ex,
    std::shared_ptr<const PhisEstablisher> phis_establisher)
    :  //_spin_realm(maths_in_physic::SpinRealm(2)),
      _spin_realm(spin_realm),
      _A(A),
      _B(B),
      _C(C),
      _D(D),
      _J(J),
      _Ez(Ez),
      _Ex(Ex),
      _phis_establisher(phis_establisher) {
  std::cout << "[DEBUG  ] [JABCDZXHamiltonianCreator] [ctor-data] ";
  std::cout << "A, B, C, D : " << _A << ", " << _B << ", " << _C << ", " << _D
            << std::endl;
  std::cout << "[DEBUG  ] [JABCDZXHamiltonianCreator] [ctor-data] ";
  std::cout << "J, Ez, Ex  : " << _J << ", " << _Ez << ", " << _Ex << std::endl;
}

void JABCDZXHamiltonianCreator::creat_hamiltonian(
    arma::mat &H, const stars_ring_basis::LocalizedBasis &basis) const {
  throw std::logic_error("not implementd!");
}

void JABCDZXHamiltonianCreator::creat_hamiltonian(
    arma::sp_mat &H, const stars_ring_basis::LocalizedBasis &basis) const {
  assert(basis.physical_system());
  const unsigned n_sites = basis.physical_system()->n_sites();
  const unsigned n_elems = basis.size();
  H = arma::sp_mat(n_elems, n_elems);
  for (unsigned ket_idx = 0; ket_idx < n_elems; ++ket_idx) {
    std::shared_ptr<stars_ring_basis::LocalizedElement> ket =
        basis.vec_index()[ket_idx];
    const stars_ring_basis::RawState ket_spin = ket->raw_state();
    const RawPhis ket_phis = _phis_establisher->establish(ket_spin);
    {
      double diag_energy = 0.0;
      for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
        const double spin_op = _A - _B * _spin_realm.S_z_diag(ket_spin[i]) *
                                        _spin_realm.S_z_diag(ket_spin[j]);
        const double orbit_op =
            psi1_orbitij_psi2(ket_phis, i, j, ket_phis, _C, _D);
        diag_energy += _J * spin_op * orbit_op;
      }
      for (unsigned i = 0; i < n_sites; i++) {
        diag_energy += -_Ez * phi1_tauc_psi2(ket_phis[i], ket_phis[i]);
        diag_energy += -_Ex * phi1_tauanty_psi2(ket_phis[i], ket_phis[i]);
      }
      H(ket_idx, ket_idx) = diag_energy;
    }
    {
      for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
        auto temp = stars_ring_basis::operations::crcr(
            i, j, ket->raw_state(), _spin_realm.multiplicity);
        if (!temp) continue;
        const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
        if (!bra_idx)
          continue;  // nie wszedl do bazy, bo baza nie musi byc zupelna!
        std::shared_ptr<stars_ring_basis::LocalizedElement> bra =
            basis.vec_index()[*bra_idx];
        const stars_ring_basis::RawState bra_spin = bra->raw_state();
        const RawPhis bra_phis = _phis_establisher->establish(bra_spin);
        const double spin_op = _B * _spin_realm.S_pm_diag(ket_spin[i]) *
                               _spin_realm.S_pm_diag(ket_spin[j]) / 2.0;
        const double orbit_op =
            psi1_orbitij_psi2(bra_phis, i, j, ket_phis, _C, _D);
        H(*bra_idx, ket_idx) = _J * spin_op * orbit_op;
      }
    }
    {
      for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
        auto temp = stars_ring_basis::operations::anan(i, j, ket->raw_state());
        if (!temp) continue;
        const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
        if (!bra_idx)
          continue;  // nie wszedl do bazy, bo baza nie musi byc zupelna!
        std::shared_ptr<stars_ring_basis::LocalizedElement> bra =
            basis.vec_index()[*bra_idx];
        const stars_ring_basis::RawState bra_spin = bra->raw_state();
        const RawPhis bra_phis = _phis_establisher->establish(bra_spin);
        const double spin_op = _B * _spin_realm.S_pm_diag(ket_spin[i] - 1) *
                               _spin_realm.S_pm_diag(ket_spin[j] - 1) / 2.0;
        const double orbit_op =
            psi1_orbitij_psi2(bra_phis, i, j, ket_phis, _C, _D);
        H(*bra_idx, ket_idx) = _J * spin_op * orbit_op;
      }
    }
  }
  // arma::mat temp_H(H);
  // temp_H(arma::span(0, 8), arma::span(0, 8)).print("H");
}

}  // namespace stars_ring_numerical
