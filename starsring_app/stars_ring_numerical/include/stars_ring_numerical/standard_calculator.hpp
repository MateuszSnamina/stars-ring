#ifndef STANDARD_JOB_CALCULATOR_HPP
#define STANDARD_JOB_CALCULATOR_HPP

#include<memory>
#include<armadillo>

#include<stars_ring_basis/basis.hpp>
#include<stars_ring_core/physical_system.hpp>
#include<stars_ring_basis/basis_toolbox.hpp>
#include<stars_ring_basis/raw_state_coupled_elements_generator.hpp>

namespace stars_ring_numerical {

    class HamiltonianCreator {
    public:
        virtual void creat_hamiltonian(
                arma::mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const = 0;
        virtual void creat_hamiltonian(
                arma::sp_mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const = 0;
        virtual ~HamiltonianCreator() = default;
    };

    class StandardCalculator {
    public:
        StandardCalculator(
                std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
                const std::vector<stars_ring_basis::RawState> init_states,
                std::shared_ptr<stars_ring_basis::CoupledRawStatesGenerator> coupled_raw_states_generator,
                unsigned perturbation_level,
                std::shared_ptr<HamiltonianCreator> hamiltonian_creator);
        // strategy:
        std::shared_ptr<HamiltonianCreator> hamiltonian_creator() const;
        void hamiltonian_creator(std::shared_ptr<HamiltonianCreator>);
        // actions:
        void do_calculations_sparse_diagonalize();
        void do_calculations_sparse_create_hamiltonian();
        void do_calculations_sparse(unsigned requested_n_states_to_calculate);
        void do_calculations_dense_create_af_hamiltonian();
        void do_calculations_dense_diagonalize();
        void do_calculations_dense();
        // accesors valid after construction:
        const stars_ring_basis::BasisBox & basis_box() const;
        // accesors valid after calculations:
        const arma::sp_mat & hamiltonian_sparse() const;
        const arma::mat & hamiltonian_dense() const;
        unsigned n_calculated_states() const;
        const arma::vec & energies() const;
        const arma::mat & states() const;
        std::tuple<const arma::vec&, const arma::mat&> energies_states() const;
        double energies(unsigned state_idx) const;
        const arma::vec states(unsigned state_idx) const;
        std::tuple<double, const arma::vec&> energies_states(unsigned state_idx) const;
    private:
        const stars_ring_basis::BasisBox _basis_box;
        std::shared_ptr<HamiltonianCreator> _hamiltonian_creator;
        arma::sp_mat _hamiltonian_sparse;
        arma::mat _hamiltonian_dense;
        unsigned _n_calculated_states; // used in sparse calculations;
        arma::vec _energies;
        arma::mat _states;
    };

}

#endif