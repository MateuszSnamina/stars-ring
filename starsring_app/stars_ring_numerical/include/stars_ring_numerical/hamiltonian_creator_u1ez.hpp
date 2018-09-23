#ifndef HAMILTONIAN_CREATOR_U1EZ_HPP
#define HAMILTONIAN_CREATOR_U1EZ_HPP
/*
#include<armadillo>
#include<memory>

#include<maths_in_physic/spin_realm.hpp>
#include<stars_ring_basis/basis.hpp>
#include<stars_ring_numerical/phis_establisher.hpp>
#include<stars_ring_numerical/standard_calculator.hpp>

namespace stars_ring_numerical {

    class U1EzHamiltonianCreator : public HamiltonianCreator {
    public:
        U1EzHamiltonianCreator(
                double E_U1, double Ez,
                //maths_in_physic::SpinRealm spin_realm,
                std::shared_ptr<const PhisEstablisher> phis_establisher);
        virtual void creat_hamiltonian(
                arma::mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const override;
        virtual void creat_hamiltonian(
                arma::sp_mat & H,
                const stars_ring_basis::LocalizedBasis& basis) const override;
    private:
        double _E_U1, _Ez;
        const maths_in_physic::SpinRealm spin_realm;
        std::shared_ptr<const PhisEstablisher> _phis_establisher;
    };

}
*/

#endif
