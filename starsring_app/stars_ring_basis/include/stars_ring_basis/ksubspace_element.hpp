#ifndef KSUBSPACE_ELEMENT_HPP
#define KSUBSPACE_ELEMENT_HPP

#include <iostream>
#include <memory>

#include <stars_ring_basis/translation_unique_element.hpp>
#include <stars_ring_core/physical_system.hpp>

namespace stars_ring_basis {

class KsubspaceElement : public stars_ring_core::SettledInPhysicalSystem {
 public:
  // Physical properties getters:
  const std::shared_ptr<TranslationUniqueElement> translation_unique_element()
      const;
  unsigned nk() const;
  double norm() const;
  // Factory:
  friend std::unique_ptr<KsubspaceElement> make_ksubspace_element(
      std::shared_ptr<TranslationUniqueElement> translation_unique_element,
      unsigned nk);
  // Member function that enables
  // a class instance to be stored in VecMap
  typedef const RawState& KeyT;
  const KeyT key() const;

 private:
  KsubspaceElement(
      std::shared_ptr<TranslationUniqueElement> translation_unique_element,
      unsigned nk);
  const std::shared_ptr<TranslationUniqueElement> _translation_unique_element;
  const unsigned _nk;
};

std::ostream& operator<<(std::ostream& s, const KsubspaceElement& el);

// Inline functions implementations:

inline const std::shared_ptr<TranslationUniqueElement>
KsubspaceElement::translation_unique_element() const {
  return _translation_unique_element;
}

inline unsigned KsubspaceElement::nk() const { return _nk; }

inline typename KsubspaceElement::KeyT KsubspaceElement::key() const {
  return translation_unique_element()->key();
}

inline double KsubspaceElement::norm() const {
  return translation_unique_element()->norm(_nk);
}

inline std::unique_ptr<KsubspaceElement> make_ksubspace_element(
    std::shared_ptr<TranslationUniqueElement> translation_unique_element,
    unsigned nk) {
  std::unique_ptr<KsubspaceElement> result(
      new KsubspaceElement(translation_unique_element, nk));
  return result;
}

}  // namespace stars_ring_basis

#endif