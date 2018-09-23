#ifndef BASIS_TYPEDEFS_HPP
#define BASIS_TYPEDEFS_HPP

#include <stars_ring_basis/basis.hpp>

namespace stars_ring_basis {

typedef Basis<LocalizedElement> LocalizedBasis;
typedef Basis<TranslationUniqueElement> TranslationUniqueBasis;
typedef Basis<KsubspaceElement> KsubspaceBasis;

}  // namespace stars_ring_basis

#endif