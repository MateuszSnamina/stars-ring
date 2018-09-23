#ifndef PHYSICAL_SYSTEM_HPP
#define PHYSICAL_SYSTEM_HPP

#include<cassert>
#include<memory>
#include<complex>
#include<iostream>

namespace stars_ring_core {

    class PhysicalSystem {
    public:
        unsigned n_sites() const;
        unsigned n_cells() const;
        unsigned n_sites_in_cell() const;
        unsigned max_n_stars() const;
        double k_unit1() const;
        double k_unit2() const;
        std::complex<double> chi() const;
        class Builder;
    private:
        PhysicalSystem(unsigned n_cells, unsigned max_n_stars);
        const unsigned _n_sites;
        const unsigned _n_cells;
        const unsigned _n_sites_in_cell;
        const unsigned _max_n_stars;
        const double _k_unit1;
        const double _k_unit2;
        const std::complex<double> _chi;
    };

    class PhysicalSystem::Builder {
    public:
        Builder();
        Builder& n_cells(unsigned n_cells);
        Builder& max_n_stars(unsigned max_n_stars);
        std::unique_ptr<PhysicalSystem> build();
    private:
        unsigned _n_cells;
        unsigned _max_n_stars;
    };

    std::ostream& operator<<(std::ostream&, const PhysicalSystem&);

    class SettledInPhysicalSystem {
    public:
        SettledInPhysicalSystem(std::shared_ptr<PhysicalSystem> physical_system);
        std::shared_ptr<PhysicalSystem> physical_system() const;
    private:
        const std::shared_ptr<PhysicalSystem> _physical_system;
    };


    // Inlined functions:

    inline
    unsigned PhysicalSystem::n_sites() const {
        return _n_sites;
    }

    inline
    unsigned PhysicalSystem::n_cells() const {
        return _n_cells;
    }

    inline
    unsigned PhysicalSystem::n_sites_in_cell() const {
        return _n_sites_in_cell;
    }

    inline
    unsigned PhysicalSystem::max_n_stars() const {
        return _max_n_stars;
    }

    inline
    double PhysicalSystem::k_unit1() const {
        return _k_unit1;
    }

    inline
    double PhysicalSystem::k_unit2() const {
        return _k_unit2;
    }

    inline
    std::complex<double> PhysicalSystem::chi() const {
        return _chi;
    }

    inline
    std::shared_ptr<PhysicalSystem>
    SettledInPhysicalSystem::physical_system() const {
        return _physical_system;
    }

}

#endif