#ifndef TRANSLATION_UNIQUE_ELEMENT_HPP
#define TRANSLATION_UNIQUE_ELEMENT_HPP

#include<memory>
#include<vector>
#include<iostream>

#include<stars_ring_core/physical_system.hpp>
#include<stars_ring_basis/raw_state.hpp>

namespace stars_ring_basis {

    class TranslationUniqueElement : public stars_ring_core::SettledInPhysicalSystem {
    public:
        // Physical properties getters:
        const RawState & raw_representative_state() const;
        const std::vector<RawState> raw_equivalent_states() const;
        unsigned cycle_length() const;
        unsigned perturbation_order() const;
        unsigned n_stars_A() const;
        unsigned n_stars_B() const;
        unsigned n_stars() const;
        bool support_nk(unsigned nk) const;
        double norm(unsigned nk) const;
        // Factory:
        friend
        std::unique_ptr<TranslationUniqueElement>
        make_translation_unique_element(
                std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
                const RawState & state,
                unsigned perturbation_order);
        // Member function that enables
        // a class instance to be stored in VecMap
        typedef const RawState& KeyT;
        const KeyT key() const;
    private:
        TranslationUniqueElement(
                std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
                const std::vector<unsigned> & state,
                unsigned perturbation_order);
        const RawState _raw_representative_state;
        //const std::vector<RawState> _raw_equivalent_states;
        const unsigned _cycle_length;
        const unsigned _perturbation_order;
        const unsigned _n_stars_A;
        const unsigned _n_stars_B;
        const unsigned _n_stars;
    };

    std::unique_ptr<TranslationUniqueElement> make_translation_unique_element(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            const std::vector<unsigned> & state,
            unsigned perturbation_order);

    std::ostream& operator<<(
            std::ostream& s,
            const TranslationUniqueElement& el);

    // Inline functions implementations:

    inline
    const RawState &
    TranslationUniqueElement::raw_representative_state() const {
        return _raw_representative_state;
    }

    inline
    unsigned
    TranslationUniqueElement::cycle_length() const {
        return _cycle_length;
    }

    inline
    unsigned
    TranslationUniqueElement::perturbation_order() const {
        return _perturbation_order;
    }

    inline
    typename TranslationUniqueElement::KeyT
    TranslationUniqueElement::key() const {
        return _raw_representative_state;
    }

    inline
    unsigned
    TranslationUniqueElement::n_stars_A() const {
        return _n_stars_A;
    }

    inline
    unsigned
    TranslationUniqueElement::n_stars_B() const {
        return _n_stars_B;
    }

    inline
    unsigned
    TranslationUniqueElement::n_stars() const {
        return _n_stars;
    }

}

#endif