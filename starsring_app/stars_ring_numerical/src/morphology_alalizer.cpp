#include <stars_ring_basis/basis_linkers.hpp>
#include <stars_ring_numerical/morphology_alalizer.hpp>

namespace stars_ring_numerical {

MorphologyAlalizer::MorphologyAlalizer(
    const stars_ring_basis::BasisBox& basis_box)
    : _basis_box(basis_box),
      _linker(construct_ksubspace_basis_to_localized_basis_linker(basis_box)) {}

void MorphologyAlalizer::calculate_morphology(const arma::vec& x) {
  calculate_full_momentum_morphology(x);
  calculate_tantamount_momentum_morphology();
  calculate_main_nk();
  calculate_average_stars_morphology(x);
}

void MorphologyAlalizer::calculate_full_momentum_morphology(
    const arma::vec& x) {
  assert(x.n_elem == _basis_box.localized_basis().size());
  assert(_basis_box.localized_basis().physical_system());
  // std::cout << "[INFO   ] [PROGRESS] Program is about to calculate the state
  // full morphology." << std::endl; arma::wall_clock timer; timer.tic();
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      _basis_box.localized_basis().physical_system();
  _full_momentum_morphology = std::vector<double>(physical_system->n_cells());
  for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
    double content = 0.0;
    for (const stars_ring_basis::LinearCombination<
             std::complex<double>, unsigned>& linear_combination :
         _linker[nk]) {
      std::complex<double> dot_product = 0.0;
      for (const stars_ring_basis::LinearCombinationIngredient<
               std::complex<double>, unsigned>& linear_combination_ingredient :
           linear_combination)
        dot_product += linear_combination_ingredient.coeficient *
                       x(linear_combination_ingredient.vector);
      content += std::norm(dot_product);  // std::norm results the sum of
                                          // squares: norm(a+bi) = a**2 + b**2
    }
    _full_momentum_morphology[nk] = content;
  }
  // const double calculation_time = timer.toc();
  // std::cout << "[INFO   ] [PROGRESS] Program has calculated the state full
  // morphology"
  //       << " (in " << calculation_time << "s)." << std::endl;
  // for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
  //    std::cout << "[DEBUG] _full_momentum_morphology:"
  //            << "(nk = " << nk << "): "
  //            << " " << _full_momentum_morphology[nk] << std::endl;
  //}
}

void MorphologyAlalizer::calculate_tantamount_momentum_morphology() {
  // std::cout << "[INFO   ] [PROGRESS] Program is about to calculate the state
  // restricted morphology." << std::endl;
  unsigned full_morphology_size = _basis_box.physical_system()->n_cells();
  unsigned restricted_morphology_size =
      _basis_box.physical_system()->n_cells() / 2 + 1;
  _tantamount_momentum_morphology =
      std::vector<double>(restricted_morphology_size);
  for (unsigned i = 0; i < restricted_morphology_size; i++) {
    _tantamount_momentum_morphology[i] += _full_momentum_morphology[i];
  }
  for (unsigned i = restricted_morphology_size; i < full_morphology_size; i++) {
    assert(full_morphology_size - i < restricted_morphology_size);
    _tantamount_momentum_morphology[full_morphology_size - i] +=
        _full_momentum_morphology[i];
  }
  // std::cout << "[INFO   ] [PROGRESS] Program has calculated the state
  // restricted morphology." << std::endl; for (unsigned nk = 0; nk <
  // restricted_morphology_size; nk++) {
  //    std::cout << "[DEBUG] _tantamount_momentum_morphology, content"
  //            << "(nk = " << nk << "): "
  //            << " " << _tantamount_momentum_morphology[nk] << std::endl;
  //}
}

void MorphologyAlalizer::calculate_main_nk() {
  for (unsigned nk = 0; nk < _tantamount_momentum_morphology.size(); nk++) {
    if (_tantamount_momentum_morphology[nk] > 0.99) {
      _main_nk = boost::optional<unsigned>(nk);
      // std::cout << "[DEBUG] found main nk = " << nk << std::endl;
      return;
    }
  }
  _main_nk = boost::optional<unsigned>();
}

void MorphologyAlalizer::calculate_average_stars_morphology(
    const arma::vec& x) {
  assert(x.n_elem == _basis_box.localized_basis().size());
  _average_n_stars_A = 0;
  _average_n_stars_B = 0;
  _average_n_stars = 0;
  for (unsigned i = 0; i < _basis_box.localized_basis().size(); i++) {
    const unsigned n_stars_A =
        _basis_box.localized_basis().vec_index()[i]->n_stars_A();
    const unsigned n_stars_B =
        _basis_box.localized_basis().vec_index()[i]->n_stars_B();
    const unsigned n_stars =
        _basis_box.localized_basis().vec_index()[i]->n_stars();
    const double probability = x(i) * x(i);
    _average_n_stars_A += probability * n_stars_A;
    _average_n_stars_B += probability * n_stars_B;
    _average_n_stars += probability * n_stars;
  }
}

const std::vector<double>& MorphologyAlalizer::full_momentum_morphology()
    const {
  return _full_momentum_morphology;
}

const std::vector<double>& MorphologyAlalizer::tantamount_momentum_morphology()
    const {
  return _tantamount_momentum_morphology;
}

const boost::optional<unsigned>& MorphologyAlalizer::main_nk() const {
  return _main_nk;
}

double MorphologyAlalizer::average_n_stars_A() const {
  return _average_n_stars_A;
}

double MorphologyAlalizer::average_n_stars_B() const {
  return _average_n_stars_B;
}

double MorphologyAlalizer::average_n_stars() const { return _average_n_stars; }

// ----------------------------------

DegeneracySubspacesAnalizer::DegeneracySubspacesAnalizer(
    const arma::vec& original_energies, const arma::mat& original_beta,
    const stars_ring_basis::BasisBox& basis_box, double energy_threshold)
    : _original_energies(original_energies),
      _original_beta(original_beta),
      _basis_box(basis_box),
      _linker(construct_ksubspace_basis_to_localized_basis_linker(basis_box)),
      _energy_threshold(energy_threshold) {
  assert(_original_beta.n_rows == basis_box.localized_basis().size());
  assert(_original_beta.n_cols <= _original_beta.n_rows);
  assert(_original_beta.n_cols == _original_energies.n_elem);
}

void DegeneracySubspacesAnalizer::determine_original_states_content() {
  assert(_basis_box.localized_basis().physical_system());
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      _basis_box.localized_basis().physical_system();
  _content_original_states = std::vector<std::vector<double>>(
      _original_beta.n_cols,
      std::vector<double>(physical_system->n_cells(), arma::datum::nan));
  for (unsigned n_state = 0; n_state < _original_beta.n_cols; ++n_state) {
    arma::vec x = _original_beta.col(n_state);
    assert(x.n_elem == _basis_box.localized_basis().size());
    for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
      double content = 0.0;
      for (const stars_ring_basis::LinearCombination<
               std::complex<double>, unsigned>& linear_combination :
           _linker[nk]) {
        std::complex<double> dot_product = 0.0;
        for (const stars_ring_basis::LinearCombinationIngredient<
                 std::complex<double>, unsigned>&
                 linear_combination_ingredient : linear_combination)
          dot_product += linear_combination_ingredient.coeficient *
                         x(linear_combination_ingredient.vector);
        content +=
            std::norm(dot_product);  // std::norm (from <complex>) results the
                                     // sum of squares: norm(a+bi) = a**2 + b**2
      }
      _content_original_states[n_state][nk] = content;
    }
  }
}

void DegeneracySubspacesAnalizer::determine_subspaces_indices() {
  unsigned current_from_idx = 0;
  double energy = _original_energies(0);
  for (unsigned idx = 1; idx < _original_energies.n_rows; ++idx) {
    assert(idx > current_from_idx);
    assert(_original_energies(idx) - _original_energies(current_from_idx) >= 0);
    if (_original_energies(idx) - _original_energies(current_from_idx) >
        _energy_threshold) {
      _subspace_indices.push_back(std::make_pair(current_from_idx, idx));
      energy /= idx - current_from_idx;
      _subspace_energies.push_back(energy);
      current_from_idx = idx;
      energy = 0.0;
    }
    energy += _original_energies(idx);
  }
  _subspace_indices.push_back(
      std::make_pair(current_from_idx, _original_energies.n_rows));
  energy /= _original_energies.n_rows - current_from_idx;
  _subspace_energies.push_back(energy);
}

void DegeneracySubspacesAnalizer::determine_subspaces_content() {
  assert(_basis_box.localized_basis().physical_system());
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      _basis_box.localized_basis().physical_system();
  _content_subspaces = std::vector<std::vector<double>>(
      _subspace_indices.size(),
      std::vector<double>(physical_system->n_cells(), 0.0));
  for (unsigned n_subspace = 0; n_subspace < _subspace_indices.size();
       ++n_subspace)
    for (unsigned nk = 0; nk < physical_system->n_cells(); nk++)
      for (unsigned n_state = _subspace_indices[n_subspace].first;
           n_state < _subspace_indices[n_subspace].second; ++n_state)
        _content_subspaces[n_subspace][nk] +=
            _content_original_states[n_state][nk];
}

void DegeneracySubspacesAnalizer::determine_subspaces_validity() {
  assert(_basis_box.localized_basis().physical_system());
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      _basis_box.localized_basis().physical_system();
  _subspaces_validity = std::vector<bool>(_subspace_indices.size(), true);
  for (unsigned n_subspace = 0; n_subspace < _content_subspaces.size();
       ++n_subspace) {
    for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
      const unsigned content = std::round(_content_subspaces[n_subspace][nk]);
      if (std::abs(content - _content_subspaces[n_subspace][nk]) > 1e-6) {
        _subspaces_validity[n_subspace] = false;
        std::cout << "[WARNING] [CONTENT FILTER] problem with n_subspace, nk, "
                     "content:"
                  << n_subspace << " " << nk << " "
                  << _content_subspaces[n_subspace][nk] << std::endl;
        if (n_subspace + 1 != _content_subspaces.size()) {
          std::cerr << "[ERROR  ] [CONTENT FILTER] error while processing a "
                       "subspace that is not the last subspace.";
        }
      }
    }
  }
}

void DegeneracySubspacesAnalizer::determine_filtered_subspaces_content() {
  assert(_basis_box.localized_basis().physical_system());
  std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system =
      _basis_box.localized_basis().physical_system();
  _fileted_content_subspaces =
      std::vector<std::vector<unsigned>>(_subspace_indices.size());
  for (unsigned n_subspace = 0; n_subspace < _content_subspaces.size();
       ++n_subspace) {
    if (_subspaces_validity[n_subspace]) {
      for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
        unsigned content = std::round(_content_subspaces[n_subspace][nk]);
        while (content-- != 0) {
          _fileted_content_subspaces[n_subspace].push_back(nk);
        }
      }
    }
  }
}

// void DegeneracySubspacesAnalizer::determine_adapted_basis() {
//   // to do.
// }

unsigned DegeneracySubspacesAnalizer::n_original_states() const {
  return _original_energies.n_rows;
}

arma::vec DegeneracySubspacesAnalizer::original_state_energies() const {
  return _original_energies;
}

double DegeneracySubspacesAnalizer::original_state_energy(
    unsigned n_state) const {
  assert(n_state < _original_energies.n_rows);
  return _original_energies(n_state);
}

unsigned DegeneracySubspacesAnalizer::n_subspaces() const {
  return _subspace_indices.size();
}

std::vector<std::pair<unsigned, unsigned>>
DegeneracySubspacesAnalizer::subspace_indices() const {
  return _subspace_indices;
}

std::pair<unsigned, unsigned> DegeneracySubspacesAnalizer::subspace_indices(
    unsigned n_subspaces) const {
  assert(n_subspaces < _subspace_indices.size());
  return _subspace_indices[n_subspaces];
}

std::vector<double> DegeneracySubspacesAnalizer::subspace_energies() const {
  return _subspace_energies;
}

double DegeneracySubspacesAnalizer::subspace_energy(
    unsigned n_subspaces) const {
  assert(n_subspaces < _subspace_energies.size());
  return _subspace_energies[n_subspaces];
}

std::vector<std::vector<double>>
DegeneracySubspacesAnalizer::content_original_states() const {
  return _content_original_states;
}

std::vector<double> DegeneracySubspacesAnalizer::content_original_state(
    unsigned n_state) const {
  assert(n_state < _content_original_states.size());
  return _content_original_states[n_state];
}

std::vector<std::vector<double>>
DegeneracySubspacesAnalizer::content_subspaces() const {
  return _content_subspaces;
}

std::vector<double> DegeneracySubspacesAnalizer::content_subspace(
    unsigned n_subspaces) const {
  assert(n_subspaces < _content_subspaces.size());
  return _content_subspaces[n_subspaces];
}

std::vector<bool> DegeneracySubspacesAnalizer::subspaces_validity() const {
  return _subspaces_validity;
}

bool DegeneracySubspacesAnalizer::subspace_validity(
    unsigned n_subspaces) const {
  assert(n_subspaces < _subspaces_validity.size());
  return _subspaces_validity[n_subspaces];
}

std::vector<std::vector<unsigned>>
DegeneracySubspacesAnalizer::fileted_content_subspaces() const {
  return _fileted_content_subspaces;
}

std::vector<unsigned> DegeneracySubspacesAnalizer::fileted_content_subspace(
    unsigned n_subspaces) const {
  assert(n_subspaces < _fileted_content_subspaces.size());
  return _fileted_content_subspaces[n_subspaces];
}

// arma::mat DegeneracySubspacesAnalizer::adapted_basis() const {
//   return _adapted_basis;
// }

}  // namespace stars_ring_numerical