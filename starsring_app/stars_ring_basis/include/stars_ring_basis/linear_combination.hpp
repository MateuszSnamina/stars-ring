#ifndef LINEAR_COMBINATION_HPP
#define LINEAR_COMBINATION_HPP

namespace stars_ring_basis {

template <typename CoefT, typename VectorT>
struct LinearCombinationIngredient {
  CoefT coeficient;
  VectorT vector;
};

template <typename CoefT, typename VectorT>
using LinearCombination =
    std::vector<LinearCombinationIngredient<CoefT, VectorT>>;

}  // namespace stars_ring_basis

#endif