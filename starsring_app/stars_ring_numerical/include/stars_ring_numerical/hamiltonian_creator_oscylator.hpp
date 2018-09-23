#ifndef OSCYLATOR_HAMILTONIAN_CREATOR_HPP
#define OSCYLATOR_HAMILTONIAN_CREATOR_HPP

#include <armadillo>

#include <maths_in_physic/oscilator_realm.hpp>
#include <stars_ring_basis/basis.hpp>
#include <stars_ring_numerical/standard_calculator.hpp>

namespace stars_ring_numerical {

class OscylatorHamiltonianCreator : public HamiltonianCreator {
 public:
  OscylatorHamiltonianCreator(maths_in_physic::OscilatorRealm oscilator_realm);
  virtual void creat_hamiltonian(
      arma::mat& H,
      const stars_ring_basis::LocalizedBasis& basis) const override;
  virtual void creat_hamiltonian(
      arma::sp_mat& H,
      const stars_ring_basis::LocalizedBasis& basis) const override;

 private:
  const maths_in_physic::OscilatorRealm oscilator_realm;
};

}  // namespace stars_ring_numerical

#endif