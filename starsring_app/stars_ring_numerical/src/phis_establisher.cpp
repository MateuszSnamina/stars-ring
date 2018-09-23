#include<cassert>
#include<stars_ring_numerical/phis_establisher.hpp>

namespace stars_ring_numerical {

    TrivialPhisEstablisher::TrivialPhisEstablisher(double phi_0) :
    _phi_0(phi_0) {
    }

    RawPhis TrivialPhisEstablisher::establish(
            const stars_ring_basis::RawState& state) const {
        //for (unsigned i = 0; i < state.size(); i++)
        //assert(state[i] == 0 || state[i] == 1);
        return RawPhis(state.size(), _phi_0);
    }

    SimplePhisEstablisher::SimplePhisEstablisher(double phi_0, double delta_phi) :
    _phi_0(phi_0), _delta_phi(delta_phi) {
    }

    RawPhis SimplePhisEstablisher::establish(
            const stars_ring_basis::RawState& state) const {
        RawPhis result(state.size());
        for (unsigned i = 0; i < state.size(); i++)
            result[i] = _phi_0 + state[i] * _delta_phi;
        return result;
    }

}