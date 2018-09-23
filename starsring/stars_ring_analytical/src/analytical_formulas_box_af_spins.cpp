#include<stars_ring_analytical/analytical_formulas_box_af_spins.hpp>

namespace stars_ring_analytical {

    AnalyticalFormulasBoxAfSpins::AnalyticalFormulasBoxAfSpins(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
            unsigned multiplicity) :
    AnalyticalFormulasBox(physical_system),
    _multiplicity(multiplicity),
    _analytical_formulas_box_af_oscylators(physical_system) {
    }


    double AnalyticalFormulasBoxAfSpins::ground_state_classical_energy() const {
        const double S = (_multiplicity - 1) / 2.0;
        return  -(physical_system()->n_sites() * S * S);
    }

    double AnalyticalFormulasBoxAfSpins::ground_state_correlation_energy() const {
        const double S = (_multiplicity - 1) / 2.0;
        return S * _analytical_formulas_box_af_oscylators.ground_state_correlation_energy();
    }

    double AnalyticalFormulasBoxAfSpins::exc_state_relative_energy(unsigned nk) const {
        const double S = (_multiplicity - 1) / 2.0;
        return S * _analytical_formulas_box_af_oscylators.exc_state_relative_energy(nk);
    }
}