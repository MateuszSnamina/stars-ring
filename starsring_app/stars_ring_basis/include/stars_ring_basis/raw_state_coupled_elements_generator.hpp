#ifndef RAW_STATE_COUPLED_ELEMENTS_GENERATOR_HPP
#define RAW_STATE_COUPLED_ELEMENTS_GENERATOR_HPP

#include<vector>

#include<stars_ring_basis/raw_state.hpp>

namespace stars_ring_basis {

    class CoupledRawStatesGenerator {
    public:
        virtual std::vector<RawState> generate(const RawState & init_state) const = 0;
        virtual ~CoupledRawStatesGenerator() = default;
    };

    class CoupledRawStatesGenerator_AF : public CoupledRawStatesGenerator {
    public:
        explicit CoupledRawStatesGenerator_AF(unsigned max_n_stars);
        std::vector<RawState> generate(const RawState & init_state) const override;
    private:
        const unsigned _max_n_stars;
    };
}

#endif