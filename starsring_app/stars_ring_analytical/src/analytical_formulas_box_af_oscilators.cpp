#include<cmath>

#include<stars_ring_analytical/analytical_formulas_box_af_oscilators.hpp>

namespace stars_ring_analytical {

    AnalyticalFormulasBoxAfOscylators::AnalyticalFormulasBoxAfOscylators(
            std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system) :
    AnalyticalFormulasBox(physical_system) {
    }

    double AnalyticalFormulasBoxAfOscylators::ground_state_classical_energy() const {
        return 0.0;
    }

    double AnalyticalFormulasBoxAfOscylators::ground_state_correlation_energy() const {
        double results = -double(physical_system()->n_sites());
        for (unsigned nk = 0; nk < physical_system()->n_cells(); ++nk)
            results += exc_state_relative_energy(nk);
        return results;
    }

    double AnalyticalFormulasBoxAfOscylators::exc_state_relative_energy(unsigned nk) const {
        const double k = nk * physical_system()->k_unit1();
        return 2.0 * std::abs(std::sin(k));
    }

}