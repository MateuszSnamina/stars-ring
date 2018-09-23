#ifndef PHIS_CREATOR_HPP
#define PHIS_CREATOR_HPP

#include<vector>
#include<stars_ring_basis/raw_state.hpp>

typedef std::vector<double> RawPhis;

namespace stars_ring_numerical {

    class PhisEstablisher {
    public:
        virtual RawPhis establish(const stars_ring_basis::RawState& state) const = 0;
    };

    class TrivialPhisEstablisher : public PhisEstablisher {
    public:
        TrivialPhisEstablisher(double phi_0);
        RawPhis establish(const stars_ring_basis::RawState& state) const override;
    private:
        double _phi_0;
    };

    class SimplePhisEstablisher : public PhisEstablisher {
    public:
        SimplePhisEstablisher(double phi_0, double delta_phi);
        RawPhis establish(const stars_ring_basis::RawState& state) const override;
    private:
        double _phi_0, _delta_phi;
    };

}

#endif