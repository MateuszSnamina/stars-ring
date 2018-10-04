#ifndef MORPHOLOGY_ALALIZER_HPP
#define MORPHOLOGY_ALALIZER_HPP

#include <armadillo>
#include <boost/optional.hpp>
#include <complex>
#include <vector>

#include <stars_ring_basis/basis_toolbox.hpp>
#include <stars_ring_basis/linear_combination.hpp>

namespace stars_ring_numerical {

class MorphologyAlalizer {
 public:
  explicit MorphologyAlalizer(const stars_ring_basis::BasisBox& basis_box);
  void calculate_morphology(const arma::vec& x);
  const std::vector<double>& full_momentum_morphology() const;
  const std::vector<double>& tantamount_momentum_morphology() const;
  const boost::optional<unsigned>& main_nk() const;
  double average_n_stars_A() const;
  double average_n_stars_B() const;
  double average_n_stars() const;

 private:
  const stars_ring_basis::BasisBox _basis_box;
  std::vector<std::vector<
      stars_ring_basis::LinearCombination<std::complex<double>, unsigned>>>
      _linker;
  void calculate_full_momentum_morphology(const arma::vec& x);
  void calculate_tantamount_momentum_morphology();
  void calculate_main_nk();
  void calculate_average_stars_morphology(const arma::vec& x);
  std::vector<double> _full_momentum_morphology;
  std::vector<double> _tantamount_momentum_morphology;
  boost::optional<unsigned> _main_nk;
  double _average_n_stars_A;
  double _average_n_stars_B;
  double _average_n_stars;
};

class DegeneracySubspacesAnalizer {
 public:
  DegeneracySubspacesAnalizer(const arma::vec& original_energies,
                              const arma::mat& original_beta,
                              const stars_ring_basis::BasisBox& basis_box,
                              double energy_threshold = 1e-6);
  void determine_original_states_content();
  void determine_subspaces_indices();
  void determine_subspaces_content();
  void determine_subspaces_validity();
  void determine_filtered_subspaces_content();
  // void determine_adapted_basis();
  unsigned n_original_states() const;
  arma::vec original_state_energies() const;
  double original_state_energy(unsigned n_state) const;
  unsigned n_subspaces() const;
  std::vector<std::pair<unsigned, unsigned>> subspace_indices() const;
  std::pair<unsigned, unsigned> subspace_indices(unsigned n_subspaces) const;
  std::vector<double> subspace_energies() const;
  double subspace_energy(unsigned n_subspaces) const;
  std::vector<std::vector<double>> content_original_states() const;
  std::vector<double> content_original_state(unsigned n_state) const;
  std::vector<std::vector<double>> content_subspaces() const;
  std::vector<double> content_subspace(unsigned n_subspaces) const;
  std::vector<bool> subspaces_validity() const;
  bool subspace_validity(unsigned n_subspaces) const;
  std::vector<std::vector<unsigned>> fileted_content_subspaces() const;
  std::vector<unsigned> fileted_content_subspace(unsigned n_subspaces) const;
  // arma::mat adapted_basis() const;

 private:
  const arma::vec& _original_energies;
  const arma::mat& _original_beta;
  const stars_ring_basis::BasisBox _basis_box;
  std::vector<std::vector<
      stars_ring_basis::LinearCombination<std::complex<double>, unsigned>>>
      _linker;
  const double _energy_threshold;
  std::vector<std::pair<unsigned, unsigned>> _subspace_indices;
  std::vector<std::vector<double>>
      _content_original_states;  // first idx: n_state, second_idx: n_k
  std::vector<std::vector<double>>
      _content_subspaces;  // first idx: n_subspace, second_idx: n_k
  std::vector<bool> _subspaces_validity;  // first idx: n_subspace
  std::vector<std::vector<unsigned>>
      _fileted_content_subspaces;  // first idx: n_subspace, second_idx: the
                                   // n_state within the subspace. Values: nk.
  std::vector<double> _subspace_energies;
  //   arma::mat _adapted_basis;
};

}  // namespace stars_ring_numerical

#endif