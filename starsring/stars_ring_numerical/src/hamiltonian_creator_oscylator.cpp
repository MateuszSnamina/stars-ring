#include<stdexcept>

#include<stars_ring_basis/raw_state_operations.hpp>
#include<stars_ring_numerical/hamiltonian_creator_oscylator.hpp>

namespace stars_ring_numerical {

    OscylatorHamiltonianCreator::OscylatorHamiltonianCreator(maths_in_physic::OscilatorRealm oscilator_realm)
    : oscilator_realm(oscilator_realm) {
    }

    void OscylatorHamiltonianCreator::creat_hamiltonian(arma::mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        throw std::logic_error("not implemented");
    }

    void OscylatorHamiltonianCreator::creat_hamiltonian(arma::sp_mat & H, const stars_ring_basis::LocalizedBasis& basis) const {
        assert(basis.physical_system());
        const unsigned n_sites = basis.physical_system()->n_sites();
        const unsigned n_elems = basis.size();
        H = arma::sp_mat(n_elems, n_elems);
        //unsigned debug_licznik = 0; //debug!!!
        for (unsigned ket_idx = 0; ket_idx < n_elems; ++ket_idx) {
            std::shared_ptr<stars_ring_basis::LocalizedElement> ket = basis.vec_index()[ket_idx];
            {
                double diag_energy = 0.0;
                for (unsigned i = 0; i < n_sites; i++) {
                    diag_energy += 2 * oscilator_realm.n_diag(ket->raw_state()[i]);
                }
                H(ket_idx, ket_idx) = diag_energy;
                //debug_licznik++; //debug!!!
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    auto temp = stars_ring_basis::operations::crcr(i, j, ket->raw_state(), oscilator_realm.n_max_stars);
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    H(*bra_idx, ket_idx) = oscilator_realm.cr_an_diag(ket->raw_state()[i]) * oscilator_realm.cr_an_diag(ket->raw_state()[j]);
                    //debug_licznik++; //debug!!!                    
                }
            }
            {
                for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
                    auto temp = stars_ring_basis::operations::anan(i, j, ket->raw_state());
                    if (!temp) continue;
                    const auto bra_idx = basis.find_element_and_get_its_ra_index(*temp);
                    if (!bra_idx) continue; // nie wszedl do bazy, bo baza nie musi byc zupelna!
                    H(*bra_idx, ket_idx) = oscilator_realm.cr_an_diag(ket->raw_state()[i] - 1) * oscilator_realm.cr_an_diag(ket->raw_state()[j] - 1);
                    //debug_licznik++; //debug!!!
                }
            }
        }
        //std::cout << "debug_licznik in OscylatorHamiltonianCreator::creat_hamiltonian: " << debug_licznik << std::endl; //debug
    }

}