#ifndef STANDARD_CALCULATOR_HPP
#define STANDARD_CALCULATOR_HPP

#include<memory>
#include<vector>

#include<stars_ring_core/physical_system.hpp>

namespace stars_ring_analytical {

    class AnalyticalFormulasBox : public stars_ring_core::SettledInPhysicalSystem {
    public:
        AnalyticalFormulasBox(std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system);
        virtual double ground_state_classical_energy() const = 0;
        virtual double ground_state_correlation_energy() const = 0;
        virtual double ground_state_energy() const;
        virtual double exc_state_relative_energy(unsigned nk) const = 0;
        virtual ~AnalyticalFormulasBox() = default;
    };

    class StandardCalculator {
    public:
        StandardCalculator(std::shared_ptr<AnalyticalFormulasBox> formulas_box);
        // strategy:
        void formulas_box(std::shared_ptr<AnalyticalFormulasBox> formulas_box);
        std::shared_ptr<AnalyticalFormulasBox> formulas_box() const;
        // action:
        void calculate();
        // result accesors:
        double ground_state_classical_energy() const;
        double ground_state_correlation_energy() const;
        double ground_state_energy() const;
        unsigned n_exc_states() const;
        double exc_state_relative_energy(unsigned idx) const;
        double exc_state_absolute_energy(unsigned idx) const;
    private:
        std::shared_ptr<AnalyticalFormulasBox> _formulas_box;
        double _ground_state_classical_energy;
        double _ground_state_correlation_energy;
        double _ground_state_energy;
        std::vector<double> _exc_state_relative_energy;
    };

}

#endif

