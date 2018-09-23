#include<stars_ring_basis/raw_state_coupled_elements_generator.hpp>
#include<stars_ring_basis/raw_state_operations.hpp>

namespace stars_ring_basis {

    CoupledRawStatesGenerator_AF::CoupledRawStatesGenerator_AF(unsigned max_n_stars) :
    _max_n_stars(max_n_stars) {
    }

    std::vector<RawState> CoupledRawStatesGenerator_AF::generate(
            const RawState & init_state) const {
        std::vector<RawState> results;
        const unsigned n_sites = init_state.size();
        for (unsigned i = 0, j = 1; i < n_sites; i++, j = (j + 1) % n_sites) {
            const auto init_state_crcr = operations::crcr(i, j, init_state, _max_n_stars);
            if (init_state_crcr) results.push_back(*init_state_crcr);
            const auto init_state_anan = operations::anan(i, j, init_state);
            if (init_state_anan) results.push_back(*init_state_anan);
        }
        return results;
    }

}