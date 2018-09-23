#ifndef PERTURBATION_BASIS_MAKER_HPP
#define PERTURBATION_BASIS_MAKER_HPP

#include<memory>
#include<vector>

#include<stars_ring_basis/basis.hpp>

namespace stars_ring_basis {

    template<typename Element>
    class CoupledElementsGenerator {
    public:
        virtual std::vector<std::shared_ptr<Element>> generate(std::shared_ptr<Element>) const = 0;
        virtual ~CoupledElementsGenerator() = default;
    };

    template<typename Element>
    class PerturbationBasisMaker {
    public:
        PerturbationBasisMaker(
                std::shared_ptr<CoupledElementsGenerator<Element>> coupled_elements_generator);
        //strategy:
        std::shared_ptr<CoupledElementsGenerator<Element>> coupled_elements_generator() const;
        void coupled_elements_generator(
                std::shared_ptr<CoupledElementsGenerator<Element>> coupled_elements_generator);
        // action:
        void populate(Basis<Element>& basis, unsigned max_order);
    private:
        std::shared_ptr<CoupledElementsGenerator<Element>> _coupled_elements_generator;
    };

    // *************************************************************************
    // ********  Member functions definitions     ******************************
    // *************************************************************************

    template<typename Element>
    PerturbationBasisMaker<Element>::PerturbationBasisMaker(
            std::shared_ptr<CoupledElementsGenerator<Element>> coupled_elements_generator) :
    _coupled_elements_generator(coupled_elements_generator) {
    }

    template<typename Element>
    std::shared_ptr<CoupledElementsGenerator<Element>> PerturbationBasisMaker<Element>::coupled_elements_generator() const {
        return _coupled_elements_generator;
    }

    template<typename Element>
    void PerturbationBasisMaker<Element>::coupled_elements_generator(std::shared_ptr<CoupledElementsGenerator<Element>> coupled_elements_generator) {
        _coupled_elements_generator = coupled_elements_generator;
    }

    template<typename Element>
    void PerturbationBasisMaker<Element>::populate(Basis<Element>& basis, unsigned max_order) {
        unsigned last_chunk_size = basis.size();
        for (unsigned order = 0; order < max_order; order++) {
            const unsigned old_size = basis.size();
            assert(last_chunk_size <= old_size);
            for (unsigned idx = old_size - last_chunk_size; idx < old_size; idx++) {
                const std::shared_ptr<Element> el = basis.vec_index()[idx];
                const std::vector<std::shared_ptr < Element>> new_elements = _coupled_elements_generator->generate(el);
                for (std::shared_ptr<Element> new_element : new_elements)
                    basis.add_element(new_element);
            }
            const unsigned new_size = basis.size();
            last_chunk_size = new_size - old_size;
        }
    }

}
#endif