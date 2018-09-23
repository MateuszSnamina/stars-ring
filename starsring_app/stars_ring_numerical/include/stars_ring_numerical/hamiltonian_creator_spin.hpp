#ifndef SPIN_HAMILTONIAN_CREATOR_HPP
#define SPIN_HAMILTONIAN_CREATOR_HPP

#include<armadillo>

#include<stars_ring_basis/basis.hpp>
#include<maths_in_physic/spin_realm.hpp>
#include<stars_ring_numerical/standard_calculator.hpp>

namespace stars_ring_numerical {

    class SpinHamiltonianCreator : public HamiltonianCreator {
    public:
        SpinHamiltonianCreator(maths_in_physic::SpinRealm spin_realm);
        virtual void creat_hamiltonian(
                arma::mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const override;
        virtual void creat_hamiltonian(
                arma::sp_mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const override;
    private:
        const maths_in_physic::SpinRealm _spin_realm;
    };

}

#endif