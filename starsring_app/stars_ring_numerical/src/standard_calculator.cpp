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
  assert(_hamiltonian_sparse.n_cols == _hamiltonian_sparse.n_rows);
  assert(_n_calculated_states < _hamiltonian_sparse.n_cols);
  // std::cout << "XXXX " << _energies.n_elem << " " << _n_calculated_states
  //          << std::endl;
  if (!arma::eigs_sym(_energies, _states, _hamiltonian_sparse,
                      _n_calculated_states, "sa")) {
    if (_hamiltonian_sparse.n_rows < 100) {
      std::cout << "[INFO   ] [SPARSE] Armadillo failed to diagonalize the "
                   "hamiltonian!"
                << std::endl;
      std::cout << "[INFO   ] [SPARSE] [FALL-BACK DENSE] Armadillo failed to "
                   "diagonalize the "
                   "hamiltonian as a dense matrix!"
                << std::endl;
      if (!arma::eig_sym(_energies, _states, arma::mat(_hamiltonian_sparse))) {
        std::cerr << "[ERROR  ] [SPARSE] [FALL-BACK DENSE] Armadillo failed to "
                     "diagonalize the hamiltonian!"
                  << std::endl;
        std::cerr << "[ERROR  ] [SPARSE] [FALL-BACK DENSE] Program termination "
                     "(with exit code 20)."
                  << std::endl;
        exit(20);
      }
    } else {
      std::cerr << "[ERROR  ] [SPARSE] Armadillo failed to diagonalize the "
                   "hamiltonian!"
                << std::endl;
      std::cerr << "[ERROR  ] [SPARSE] Program termination (with exit code 20)."
                << std::endl;
      exit(20);
    }
  }
  if (_energies.n_elem < _n_calculated_states) {
    std::cout
        << "[INFO   ] [SPARSE] Armadillo failed to diagonalize the hamiltonian "
           "but not reporter an error!"
        << std::endl;
    std::cout << "[INFO   ] [SPARSE] [SECOND-TRY] The program is about to try "
                 "diagonalize the "
                 "hamiltonian with grater n_calculated_states requested."
              << std::endl;
    const unsigned n_calculated_states_second_try = std::min(
        _n_calculated_states + 3, _basis_box.localized_basis().size() - 1);
    std::cout << "[INFO   ] [SECOND-TRY] n_calculated_states, "
                 "n_calculated_states_second_try: "
              << _n_calculated_states << ", " << n_calculated_states_second_try
              << "." << std::endl;
    if (!arma::eigs_sym(_energies, _states, _hamiltonian_sparse,
                        n_calculated_states_second_try, "sa")) {
      std::cerr << "[ERROR  ] [SECOND-TRY] Armadillo failed to diagonalize the "
                   "hamiltonian!"
                << std::endl;
      std::cerr
          << "[ERROR  ] [SECOND-TRY] Program termination (with exit code 20)."
          << std::endl;
      exit(20);
    }
  }
  assert(_energies.n_elem == _states.n_cols);
  assert(_energies.n_elem >= _n_calculated_states);
  assert(_states.n_cols >= _n_calculated_states);
  assert(_states.n_rows == _hamiltonian_sparse.n_rows);
  _energies = _energies.rows(arma::span(0, _n_calculated_states - 1));
  _states = _states.cols(arma::span(0, _n_calculated_states - 1));
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
  std::cout << "[INFO   ] [PROGRESS] [DENSE] Program has diagonalized the "
               "hamiltonian"
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