#ifndef BASIS_HPP
#define BASIS_HPP

#include <boost/optional.hpp>
#include <cassert>
#include <memory>

#include <stars_ring_basis/ksubspace_element.hpp>
#include <stars_ring_basis/localized_element.hpp>
#include <stars_ring_basis/translation_unique_element.hpp>
#include <stars_ring_core/physical_system.hpp>
#include <utility_kit/vec_map.hpp>

namespace stars_ring_basis {

template <typename Element>
class Basis : public stars_ring_core::SettledInPhysicalSystem {
 private:
  typedef utility::VecMap<Element> VecMapT;
  VecMapT _vecmap;

 public:
  explicit Basis(
      std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system);
  void add_element(typename VecMapT::ElementPtrT c);
  const typename VecMapT::VecIndex& vec_index() const;
  const typename VecMapT::MapIndex& map_index() const;
  boost::optional<unsigned> find_element_and_get_its_ra_index(
      typename VecMapT::KeyT v) const;
  unsigned size() const;
};

// *************************************************************************
// ********  Member functions definitions     ******************************
// *************************************************************************

template <typename Element>
Basis<Element>::Basis(
    std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system)
    : stars_ring_core::SettledInPhysicalSystem(physical_system) {
  assert(physical_system);
}

template <typename Element>
void Basis<Element>::add_element(
    typename Basis<Element>::VecMapT::ElementPtrT c) {
  assert(c);
  assert(c->physical_system());
  assert(c->physical_system() == physical_system());
  _vecmap.add_element(c);
}

template <typename Element>
const typename Basis<Element>::VecMapT::VecIndex& Basis<Element>::vec_index()
    const {
  return _vecmap.vec_index();
}

template <typename Element>
const typename Basis<Element>::VecMapT::MapIndex& Basis<Element>::map_index()
    const {
  return _vecmap.map_index();
}

template <typename Element>
boost::optional<unsigned> Basis<Element>::find_element_and_get_its_ra_index(
    typename VecMapT::KeyT v) const {
  return _vecmap.find_element_and_get_its_ra_index(v);
}

template <typename Element>
unsigned Basis<Element>::size() const {
  return _vecmap.size();
}
}  // namespace stars_ring_basis

#endif