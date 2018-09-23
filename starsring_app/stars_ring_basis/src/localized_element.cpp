#include<iomanip>

#include<stars_ring_basis/raw_state_operations.hpp>
#include<stars_ring_basis/localized_element.hpp>

namespace stars_ring_basis {

    LocalizedElement::LocalizedElement(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            const RawState & state,
            unsigned perturbation_order) :
    SettledInPhysicalSystem(physical_system),
    _raw_state(state),
    _perturbation_order(perturbation_order),
    _n_stars_A(operations::n_stars_A(_raw_state)),
    _n_stars_B(operations::n_stars_B(_raw_state)),
    _n_stars(operations::n_stars(_raw_state)) {
        assert(physical_system);
        assert(state.size() == physical_system->n_sites());
    }

    std::unique_ptr<LocalizedElement> make_localized_element(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            const RawState & state,
            unsigned perturbation_order) {
        assert(physical_system);
        assert(state.size() == physical_system->n_sites());
        LocalizedElement * const raw = new LocalizedElement(physical_system, state, perturbation_order);
        return std::unique_ptr<LocalizedElement>(raw);
    }

    std::ostream& operator<<(
            std::ostream& s,
            const LocalizedElement& el) {
        auto f = s.flags();
        s << el.raw_state();
        s << " ";
        s << "〘" << std::setw(2) << el.perturbation_order() << "⚡〙";
        s.flags(f);
        return s;
    }

}