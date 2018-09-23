#include <armadillo>

#include <stars_ring_basis/build_perturbation_localized_basis.hpp>

namespace stars_ring_basis {

LocalizedBasis build_localized_basis(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    std::vector<RawState> init_states,
    std::shared_ptr<CoupledElementsGenerator<LocalizedElement>>
        coupled_localized_elements_generator,
    unsigned perturbation_level) {
  std::cout
      << "[INFO   ] [PROGRESS] Program is about to construct localized-basis."
      << std::endl;
  arma::wall_clock timer;
  timer.tic();
  // -------------------------------
  assert(physical_system);
  assert(coupled_localized_elements_generator);
  // -------------------------------
  LocalizedBasis localized_basis(physical_system);
  for (const RawState& init_state : init_states) {
    for (unsigned i = 0; i < physical_system->n_cells(); i++) {
      RawState temp = init_state;
      std::rotate(temp.begin(), temp.begin() + 2 * i, temp.end());
      localized_basis.add_element(
          make_localized_element(physical_system, temp, 0));
    }
  }
  PerturbationBasisMaker<LocalizedElement> perturbation_basis_generator(
      coupled_localized_elements_generator);
  perturbation_basis_generator.populate(localized_basis, perturbation_level);
  // -------------------------------
  const double calculations_time = timer.toc();
  std::cout << "[INFO   ] [PROGRESS] Program has constructed localized-basis"
            << " (in " << calculations_time << "s)." << std::endl;
  return localized_basis;
}

}  // namespace stars_ring_basis
