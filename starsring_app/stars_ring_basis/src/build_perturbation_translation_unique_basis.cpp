#include <armadillo>

#include <stars_ring_basis/build_perturbation_translation_unique_basis.hpp>

namespace stars_ring_basis {

TranslationUniqueBasis build_translation_unique_basis(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    std::vector<RawState> init_states,
    std::shared_ptr<CoupledElementsGenerator<TranslationUniqueElement>>
        coupled_translation_unique_elements_generator,
    unsigned perturbation_level) {
  std::cout << "[INFO   ] [PROGRESS] Program is about to construct "
               "translation-unique-basis."
            << std::endl;
  arma::wall_clock timer;
  timer.tic();
  // -------------------------------
  assert(physical_system);
  assert(coupled_translation_unique_elements_generator);
  // -------------------------------
  assert(physical_system);
  assert(coupled_translation_unique_elements_generator);
  // -------------------------------
  TranslationUniqueBasis translation_unique_basis(physical_system);
  for (const RawState& init_state : init_states) {
    translation_unique_basis.add_element(
        make_translation_unique_element(physical_system, init_state, 0));
  }
  PerturbationBasisMaker<TranslationUniqueElement> perturbation_basis_generator(
      coupled_translation_unique_elements_generator);
  perturbation_basis_generator.populate(translation_unique_basis,
                                        perturbation_level);
  // -------------------------------
  const double calculations_time = timer.toc();
  std::cout
      << "[INFO   ] [PROGRESS] Program has constructed translation-unique-basis"
      << " (in " << calculations_time << "s)." << std::endl;
  //---------------------
  return translation_unique_basis;
}

}  // namespace stars_ring_basis
