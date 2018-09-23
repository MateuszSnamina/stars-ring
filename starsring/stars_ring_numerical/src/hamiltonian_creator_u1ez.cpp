#include<stdexcept>

#include<stars_ring_basis/raw_state.hpp>

#include<stars_ring_basis/raw_state_operations.hpp>
#include<stars_ring_numerical/hamiltonian_creator_jabcdz.hpp>

/*
namespace {

    double phi1_id_psi2(double theta_1, double theta_2) {
        return std::cos((theta_1 - theta_2) / 2.0);
    }

    double phi1_tauc_psi2(double theta_1, double theta_2) {
        return -std::cos((theta_1 + theta_2) / 2.0) / 2.0;
    }

    double psi1_id_psi2(const RawPhis &ket_phis_1, const RawPhis & ket_phis_2) {
        assert(ket_phis_1.size() == ket_phis_2.size());
        double result = 1.0;
        for (unsigned idx = 0; idx < ket_phis_1.size(); ++idx)
            result *= phi1_id_psi2(ket_phis_1[idx], ket_phis_2[idx]);
        return result;
    }

    double psi1_tauci_taucj_psi2(
            const RawPhis &bra_phis,
            unsigned i, unsigned j,
            const RawPhis & ket_phis) {
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

    double orbit_op(
            double alpha, double beta,
            const RawPhis &bra_phis,
            unsigned i, unsigned j,
            const RawPhis & ket_phis) {
        assert(i != j);
        assert(i < bra_phis.size());
        assert(j < bra_phis.size());
        assert(bra_phis.size() == ket_phis.size());
        double orbit_ss = psi1_id_psi2(bra_phis, ket_phis);
        double orbit_tt = psi1_tauci_taucj_psi2(bra_phis, i, j, ket_phis);
        return alpha * orbit_ss + beta * orbit_tt;
    }

    double half_P_zeta_c_P_xi_c(
            const RawPhis &bra_phis,
            unsigned i, unsigned j,
            const RawPhis & ket_phis) {
        assert(i != j);
        assert(i < bra_phis.size());
        assert(j < bra_phis.size());
        assert(bra_phis.size() == ket_phis.size());
        return orbit_op(
                1.0 / 4.0, - 1.0,
                bra_phis, i, j, ket_phis);
    }

}

namespace stars_ring_numerical {

    U1EzHamiltonianCreator::U1EzHamiltonianCreator(
            double E_U1, double Ez,
            //maths_in_physic::SpinRealm spin_realm,
            std::shared_ptr<const PhisEstablisher> phis_establisher) :
    _E_U1(E_U1), _Ez(Ez),
    spin_realm(maths_in_physic::SpinRealm(2)),
    _phis_establisher(phis_establisher) {
    }

    void U1EzHamiltonianCreator::creat_hamiltonian(arma::mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        throw std::logic_error("not implementd!");
    }

    void U1EzHamiltonianCreator::creat_hamiltonian(arma::sp_mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        assert(basis.physical_system());
        const unsigned n_sites = basis.physical_system()->n_sites();
        const unsigned n_elems = basis.size();
        H = arma::sp_mat(n_elems, n_elems);
        for (unsigned ket_idx = 0; ket_idx < n_elems; ++ket_idx) {
            std::shared_ptr<stars_ring_basis::LocalizedElement> ket = basis.vec_index()[ket_idx];
            const stars_ring_basis::RawState ket_spin = ket->raw_state();
            const RawPhis ket_phis = _phis_establisher->establish(ket_spin);
            {
                double diag_energy = 0.0;
                assert(false);
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    double spin_op = 3.0 / 4.0 - spin_realm.S_z_diag(ket_spin[i]) * spin_realm.S_z_diag(ket_spin[j]);
                    double orbit_op = half_P_zeta_c_P_xi_c(ket_phis, i, j, ket_phis);
                    diag_energy += -2 * _E_U1 * spin_op * orbit_op;
                }
                H(ket_idx, ket_idx) = diag_energy;
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    auto temp = stars_ring_basis::operations::crcr(i, j, ket->raw_state(), spin_realm.multiplicity);
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    std::shared_ptr<stars_ring_basis::LocalizedElement> bra = basis.vec_index()[*bra_idx];
                    const stars_ring_basis::RawState bra_spin = bra->raw_state();
                    const RawPhis bra_phis = _phis_establisher->establish(bra_spin);
                    double spin_op = 3.0 / 4.0 - spin_realm.S_pm_diag(ket_spin[i]) * spin_realm.S_pm_diag(ket_spin[j]) / 2.0;
                    double orbit_op = half_P_zeta_c_P_xi_c(bra_phis, i, j, ket_phis);
                    H(*bra_idx, ket_idx) = -2 * _E_U1 * spin_op * orbit_op;
                }
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    assert(false);

                    auto temp = stars_ring_basis::operations::anan(i, j, ket->raw_state());
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    H(*bra_idx, ket_idx) = spin_realm.S_pm_diag(ket_spin[i] - 1) * spin_realm.S_pm_diag(ket_spin[j] - 1) / 2.0;
                }
            }
        }
    }

}
*/