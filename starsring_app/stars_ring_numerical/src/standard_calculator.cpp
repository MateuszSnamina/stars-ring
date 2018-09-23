#include <stars_ring_numerical/standard_calculator.hpp>

namespace stars_ring_numerical {

StandardCalculator::StandardCalculator(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system,
    const std::vector<stars_ring_basis::RawState> init_states,
    std::shared_ptr<stars_ring_basis::CoupledRawStatesGenerator>
        coupled_raw_states_generator,
    unsigned perturbation_level,
    std::shared_ptr<HamiltonianCreator> hamiltonian_creator)
    : _basis_box(stars_ring_basis::BasisBox::Builder().build(
          physical_system, init_states, coupled_raw_states_generator,
          perturbation_level)),
      _hamiltonian_creator(hamiltonian_creator) {}

//-----------------------------------

std::shared_ptr<HamiltonianCreator> StandardCalculator::hamiltonian_creator()
    const {
  return _hamiltonian_creator;
}

void StandardCalculator::hamiltonian_creator(
    std::shared_ptr<HamiltonianCreator> hamiltonian_creator) {
  _hamiltonian_creator = hamiltonian_creator;
}

//-----------------------------------

void StandardCalculator::do_calculations_sparse_create_hamiltonian() {
  assert(_basis_box.localized_basis().physical_system());
  arma::wall_clock timer;
  timer.tic();
  std::cout << "[INFO   ] [PROGRESS] [SPARSE] Program is about to build the "
               "hamiltonian."
            << std::endl;
  _hamiltonian_creator->creat_hamiltonian(_hamiltonian_sparse,
                                          _basis_box.localized_basis());
  const double calculation_time = timer.toc();
  std::cout << "[INFO   ] [PROGRESS] [SPARSE] Program has build the hamiltonian"
            << " (in " << calculation_time << "s)." << std::endl;
}

void StandardCalculator::do_calculations_sparse_diagonalize() {
  std::cout << "[INFO   ] [PROGRESS] [SPARSE] Program is about to diagonalize "
               "the hamiltonian."
            << std::endl;
  arma::wall_clock timer;
  timer.tic();
  if (!arma::eigs_sym(_energies, _states, _hamiltonian_sparse,
                      _n_calculated_states, "sa")) {
    std::cerr << "[ERROR  ] Armadillo failed to diagonalize the hamiltonian!"
              << std::endl;
    std::cerr << "[ERROR  ] Program termination (with exit code 20)."
              << std::endl;
    exit(20);
  }
  const double calculation_time = timer.toc();
  std::cout << "[INFO   ] [PROGRESS] [SPARSE] Program has diagonalized the "
               "hamiltonian"
            << " (in " << calculation_time << "s)." << std::endl;
}

void StandardCalculator::do_calculations_sparse(
    unsigned requested_n_states_to_calculate) {
  assert(_basis_box.localized_basis().physical_system());
  assert(_basis_box.localized_basis().size() > 1);
  do_calculations_sparse_create_hamiltonian();
  _n_calculated_states = std::min(requested_n_states_to_calculate,
                                  _basis_box.localized_basis().size() - 1);
  do_calculations_sparse_diagonalize();
}

//-----------------------------------

void StandardCalculator::do_calculations_dense_create_af_hamiltonian() {
  assert(_basis_box.localized_basis().physical_system());
  arma::wall_clock timer;
  timer.tic();
  std::cout << "[INFO   ] [PROGRESS] [DENSE] Program is about to build the "
               "hamiltonian."
            << std::endl;
  _hamiltonian_creator->creat_hamiltonian(_hamiltonian_dense,
                                          _basis_box.localized_basis());
  const double calculation_time = timer.toc();
  std::cout << "[INFO   ] [PROGRESS] [DENSE] Program has build the hamiltonian"
            << " (in " << calculation_time << "s)." << std::endl;
}

void StandardCalculator::do_calculations_dense_diagonalize() {
  std::cout << "[INFO   ] [PROGRESS] [DENSE] Program is about to diagonalize "
               "the hamiltonian."
            << std::endl;
  arma::wall_clock timer;
  timer.tic();
  if (!arma::eig_sym(_energies, _states, _hamiltonian_dense)) {
    std::cerr << "[ERROR  ] Armadillo failed to diagonalize the hamiltonian!"
              << std::endl;
    std::cerr << "[ERROR  ] Program termination (with exit code 20)."
              << std::endl;
    exit(20);
  }
  const double calculation_time = timer.toc();
  std::cout
      << "[INFO   ] [PROGRESS] [DENSE] Program has diagonalized the hamiltonian"
      << " (in " << calculation_time << "s)." << std::endl;
}

void StandardCalculator::do_calculations_dense() {
  assert(_basis_box.localized_basis().physical_system());
  do_calculations_dense_create_af_hamiltonian();
  _n_calculated_states = _basis_box.localized_basis().size();
  do_calculations_dense_diagonalize();
}

//-----------------------------------

const stars_ring_basis::BasisBox& StandardCalculator::basis_box() const {
  return _basis_box;
}

const arma::sp_mat& StandardCalculator::hamiltonian_sparse() const {
  return _hamiltonian_sparse;
}

const arma::mat& StandardCalculator::hamiltonian_dense() const {
  return _hamiltonian_dense;
}

unsigned StandardCalculator::n_calculated_states() const {
  return _n_calculated_states;
}

const arma::vec& StandardCalculator::energies() const { return _energies; }

const arma::mat& StandardCalculator::states() const { return _states; }

std::tuple<const arma::vec&, const arma::mat&>
StandardCalculator::energies_states() const {
  return std::make_tuple(_energies, _states);
}

double StandardCalculator::energies(unsigned state_idx) const {
  assert(state_idx < _energies.n_elem);
  return _energies(state_idx);
}

const arma::vec StandardCalculator::states(unsigned state_idx) const {
  assert(state_idx < _states.n_cols);
  return _states.col(state_idx);
}

std::tuple<double, const arma::vec&> StandardCalculator::energies_states(
    unsigned state_idx) const {
  assert(state_idx < _states.n_cols);
  assert(state_idx < _energies.n_elem);
  return std::make_tuple(_energies(state_idx), _states.col(state_idx));
}

}  // namespace stars_ring_numerical