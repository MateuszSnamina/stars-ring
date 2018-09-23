#include <iomanip>

#include <stars_ring_basis/ksubspace_element.hpp>

namespace stars_ring_basis {

KsubspaceElement::KsubspaceElement(
    std::shared_ptr<TranslationUniqueElement> translation_unique_element,
    unsigned nk)
    : SettledInPhysicalSystem(translation_unique_element->physical_system()),
      _translation_unique_element(translation_unique_element),
      _nk(nk) {
  assert(translation_unique_element);
  assert(translation_unique_element->raw_representative_state().size() % 2 ==
         0);
  assert(nk <
         translation_unique_element->raw_representative_state().size() / 2);
  assert(translation_unique_element->support_nk(nk));
}

std::ostream& operator<<(std::ostream& s, const KsubspaceElement& el) {
  auto f = s.flags();
  s << *el.translation_unique_element();
  s << "(nk = " << std::setw(2) << el.nk() << ")";
  s.flags(f);
  return s;
}

}  // namespace stars_ring_basis