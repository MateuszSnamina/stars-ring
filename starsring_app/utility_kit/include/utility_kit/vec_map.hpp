#ifndef VEC_MAP_HPP
#define VEC_MAP_HPP

#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/optional.hpp>
#include <memory>

namespace utility {

/*
 * VecMap<Element> is a map-like container and is a vec-like container
 * (simultaneously) for objects of class Element.
 *
 * The class Element requirements:
 * (1) Element must define the type Element::KeyT
 * (2) Element must define the Element::KeyT Element::key() const member
 * function. The VecMap containers uses above member function as a key for the
 * map.
 *
 * In the container objects are storred as shared_ptr to the objecs.
 */

template <typename Element>
class VecMap {
 public:
  typedef typename Element::KeyT KeyT;
  typedef std::shared_ptr<Element> ElementPtrT;

 private:
  // Tags for random-access-index and search-index;
  struct Vec;
  struct Map;
  // Container type definition -- helper typedefs:
  typedef boost::multi_index::tag<Vec> VecTagDef;
  typedef boost::multi_index::tag<Map> MapTagDef;
  typedef boost::multi_index::const_mem_fun<Element, KeyT, &Element::key>
      KayExtractorDef;
  // Container type definition -- index typedefs:
  typedef boost::multi_index::random_access<VecTagDef> VecIndexDef;
  typedef boost::multi_index::ordered_unique<MapTagDef, KayExtractorDef>
      MapIndexDef;
  // Container type definition -- final container typedef:
  typedef boost::multi_index::multi_index_container<
      ElementPtrT, boost::multi_index::indexed_by<VecIndexDef, MapIndexDef>>
      Container;
  // The container:
  Container container;

 public:
  // Container type definition -- index typedefs:
  typedef typename Container::template index<Vec>::type VecIndex;
  typedef typename Container::template index<Map>::type MapIndex;
  // Indices interfaces:
  VecIndex& vec_index();
  MapIndex& map_index();
  const VecIndex& vec_index() const;
  const MapIndex& map_index() const;
  // Member functions:
  void add_element(ElementPtrT c);
  boost::optional<unsigned> find_element_and_get_its_ra_index(KeyT v) const;
  unsigned size() const;
};

template <typename Element>
struct VecMap<Element>::Vec {};

template <typename Element>
struct VecMap<Element>::Map {};

// *************************************************************************
// ********  Member functions definitions     ******************************
// *************************************************************************

template <typename Element>
typename VecMap<Element>::VecIndex& VecMap<Element>::vec_index() {
  return container.template get<Vec>();
}

template <typename Element>
typename VecMap<Element>::MapIndex& VecMap<Element>::map_index() {
  return container.template get<Map>();
}

template <typename Element>
const typename VecMap<Element>::VecIndex& VecMap<Element>::vec_index() const {
  return container.template get<Vec>();
}

template <typename Element>
const typename VecMap<Element>::MapIndex& VecMap<Element>::map_index() const {
  return container.template get<Map>();
}

template <typename Element>
void VecMap<Element>::add_element(ElementPtrT c) {
  vec_index().push_back(c);
}

template <typename Element>
boost::optional<unsigned> VecMap<Element>::find_element_and_get_its_ra_index(
    KeyT v) const {
  auto search_iter = map_index().find(v);
  if (search_iter == map_index().end()) return boost::optional<unsigned>();
  auto ra_iter = container.template project<Vec>(search_iter);
  return std::distance(vec_index().begin(), ra_iter);
}

template <typename Element>
unsigned VecMap<Element>::size() const {
  return container.size();
}

}  // namespace utility

#endif