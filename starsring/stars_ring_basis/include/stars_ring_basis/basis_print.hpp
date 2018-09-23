#ifndef BASIS_PRINT_HPP
#define BASIS_PRINT_HPP

namespace stars_ring_basis {

    template<typename Element>
    void print_basis(const Basis<Element>& basis) {
        std::cout << "<BASIS>" << std::endl;
        for (unsigned idx = 0; idx < basis.size(); idx++) {
            std::cout << "  " << std::setw(5) << idx
                    << " : " << *basis.vec_index()[idx] << std::endl;
        }
        std::cout << "</BASIS>" << std::endl;
    }

}
#endif