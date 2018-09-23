#ifndef ANALYTICAL_FORMULAS_BOX_OSCILATORS_HPP
#define ANALYTICAL_FORMULAS_BOX_OSCILATORS_HPP

#include<memory>

#include<stars_ring_core/physical_system.hpp>
#include<stars_ring_analytical/standard_calculator.hpp>

namespace stars_ring_analytical {

    class AnalyticalFormulasBoxAfOscylators : public AnalyticalFormulasBox {
    public:
        AnalyticalFormulasBoxAfOscylators(
                std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system);
        double ground_state_classical_energy() const override;
        double ground_state_correlation_energy() const override;
        double exc_state_relative_energy(unsigned nk) const override;
    };

}
#endif

