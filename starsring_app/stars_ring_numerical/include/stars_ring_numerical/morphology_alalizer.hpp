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

}  // namespace stars_ring_numerical

#endif