#include<stars_ring_numerical/morphology_alalizer.hpp>
#include<stars_ring_basis/basis_linkers.hpp>

namespace stars_ring_numerical {

    MorphologyAlalizer::MorphologyAlalizer(const stars_ring_basis::BasisBox & basis_box) :
    _basis_box(basis_box),
    _linker(construct_ksubspace_basis_to_localized_basis_linker(basis_box)) {
    }

    void MorphologyAlalizer::calculate_morphology(const arma::vec & x) {
        calculate_full_momentum_morphology(x);
        calculate_tantamount_momentum_morphology();
        calculate_main_nk();
        calculate_average_stars_morphology(x);
    }

    void MorphologyAlalizer::calculate_full_momentum_morphology(const arma::vec & x) {
        assert(x.n_elem == _basis_box.localized_basis().size());
        assert(_basis_box.localized_basis().physical_system());
        //std::cout << "[INFO   ] [PROGRESS] Program is about to calculate the state full morphology." << std::endl;
        //arma::wall_clock timer;
        //timer.tic();
        std::shared_ptr<stars_ring_core::PhysicalSystem> physical_system = _basis_box.localized_basis().physical_system();
        _full_momentum_morphology = std::vector<double>(physical_system->n_cells());
        for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
            double content = 0.0;
            for (const stars_ring_basis::LinearCombination<std::complex<double>, unsigned>& linear_combination : _linker[nk]) {
                std::complex<double> dot_product = 0.0;
                for (const stars_ring_basis::LinearCombinationIngredient<std::complex<double>, unsigned> & linear_combination_ingredient : linear_combination)
                    dot_product += linear_combination_ingredient.coeficient * x(linear_combination_ingredient.vector);
                content += std::norm(dot_product); // std::norm results the sum of squares: norm(a+bi) = a**2 + b**2
            }
            _full_momentum_morphology[nk] = content;
        }
        //const double calculation_time = timer.toc();
        //std::cout << "[INFO   ] [PROGRESS] Program has calculated the state full morphology"
        //       << " (in " << calculation_time << "s)." << std::endl;
        //for (unsigned nk = 0; nk < physical_system->n_cells(); nk++) {
        //    std::cout << "[DEBUG] _full_momentum_morphology:"
        //            << "(nk = " << nk << "): "
        //            << " " << _full_momentum_morphology[nk] << std::endl;
        //}
    }

    void MorphologyAlalizer::calculate_tantamount_momentum_morphology() {
        //std::cout << "[INFO   ] [PROGRESS] Program is about to calculate the state restricted morphology." << std::endl;
        unsigned full_morphology_size = _basis_box.physical_system()->n_cells();
        unsigned restricted_morphology_size = _basis_box.physical_system()->n_cells() / 2 + 1;
        _tantamount_momentum_morphology = std::vector<double>(restricted_morphology_size);
        for (unsigned i = 0; i < restricted_morphology_size; i++) {
            _tantamount_momentum_morphology[i] += _full_momentum_morphology[i];
        }
        for (unsigned i = restricted_morphology_size; i < full_morphology_size; i++) {
            assert(full_morphology_size - i < restricted_morphology_size);
            _tantamount_momentum_morphology[full_morphology_size - i] += _full_momentum_morphology[i];
        }
        //std::cout << "[INFO   ] [PROGRESS] Program has calculated the state restricted morphology." << std::endl;
        //for (unsigned nk = 0; nk < restricted_morphology_size; nk++) {
        //    std::cout << "[DEBUG] _tantamount_momentum_morphology, content"
        //            << "(nk = " << nk << "): "
        //            << " " << _tantamount_momentum_morphology[nk] << std::endl;
        //}
    }

    void MorphologyAlalizer::calculate_main_nk() {
        for (unsigned nk = 0; nk < _tantamount_momentum_morphology.size(); nk++) {
            if (_tantamount_momentum_morphology[nk] > 0.99999) {
                _main_nk = boost::optional<unsigned>(nk);
                //std::cout << "[DEBUG] found main nk = " << nk << std::endl;
                return;
            }
        }
        _main_nk = boost::optional<unsigned>();
    }

    void MorphologyAlalizer::calculate_average_stars_morphology(const arma::vec & x) {
        assert(x.n_elem == _basis_box.localized_basis().size());
        _average_n_stars_A = 0;
        _average_n_stars_B = 0;
        _average_n_stars = 0;
        for (unsigned i = 0; i < _basis_box.localized_basis().size(); i++) {
            const unsigned n_stars_A = _basis_box.localized_basis().vec_index()[i]->n_stars_A();
            const unsigned n_stars_B = _basis_box.localized_basis().vec_index()[i]->n_stars_B();
            const unsigned n_stars = _basis_box.localized_basis().vec_index()[i]->n_stars();
            const double probability = x(i) * x(i);
            _average_n_stars_A += probability * n_stars_A;
            _average_n_stars_B += probability * n_stars_B;
            _average_n_stars += probability * n_stars;
        }
    }

    const std::vector<double>& MorphologyAlalizer::full_momentum_morphology() const {
        return _full_momentum_morphology;
    }

    const std::vector<double>& MorphologyAlalizer::tantamount_momentum_morphology() const {
        return _tantamount_momentum_morphology;
    }

    const boost::optional<unsigned>& MorphologyAlalizer::main_nk() const {
        return _main_nk;
    }

    double MorphologyAlalizer::average_n_stars_A() const {
        return _average_n_stars_A;
    }

    double MorphologyAlalizer::average_n_stars_B() const {
        return _average_n_stars_B;
    }

    double MorphologyAlalizer::average_n_stars() const {
        return _average_n_stars;
    }

}