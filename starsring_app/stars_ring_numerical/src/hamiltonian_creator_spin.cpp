#include<stdexcept>

#include<stars_ring_basis/raw_state_operations.hpp>
#include<stars_ring_numerical/hamiltonian_creator_spin.hpp>

namespace stars_ring_numerical {

    SpinHamiltonianCreator::SpinHamiltonianCreator(maths_in_physic::SpinRealm spin_realm) :
    _spin_realm(spin_realm) {
    }

    void SpinHamiltonianCreator::creat_hamiltonian(arma::mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        const unsigned n_elems = basis.size();
        H = arma::mat(n_elems, n_elems);
        for (unsigned ket_idx = 0; ket_idx < n_elems; ++ket_idx) {
            std::shared_ptr<stars_ring_basis::LocalizedElement> ket = basis.vec_index()[ket_idx];
            for (unsigned bra_idx = 0; bra_idx < n_elems; ++bra_idx) {
                std::shared_ptr<stars_ring_basis::LocalizedElement> bra = basis.vec_index()[bra_idx];
                H(bra_idx, ket_idx) = stars_ring_basis::operations::determine_coupling(bra->raw_state(), ket->raw_state(), _spin_realm);
            }
        }
    }

    void SpinHamiltonianCreator::creat_hamiltonian(arma::sp_mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        assert(basis.physical_system());
        const unsigned n_sites = basis.physical_system()->n_sites();
        const unsigned n_elems = basis.size();
        H = arma::sp_mat(n_elems, n_elems);
        for (unsigned ket_idx = 0; ket_idx < n_elems; ++ket_idx) {
            std::shared_ptr<stars_ring_basis::LocalizedElement> ket = basis.vec_index()[ket_idx];
            {
                double diag_energy = 0.0;
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    diag_energy -= _spin_realm.S_z_diag(ket->raw_state()[i]) * _spin_realm.S_z_diag(ket->raw_state()[j]);
                }
                H(ket_idx, ket_idx) = diag_energy;
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    auto temp = stars_ring_basis::operations::crcr(i, j, ket->raw_state(), _spin_realm.multiplicity);
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    H(*bra_idx, ket_idx) = _spin_realm.S_pm_diag(ket->raw_state()[i]) * _spin_realm.S_pm_diag(ket->raw_state()[j]) / 2.0;
                }
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    auto temp = stars_ring_basis::operations::anan(i, j, ket->raw_state());
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    H(*bra_idx, ket_idx) = _spin_realm.S_pm_diag(ket->raw_state()[i] - 1) * _spin_realm.S_pm_diag(ket->raw_state()[j] - 1) / 2.0;
                }
            }
        }
    }

}