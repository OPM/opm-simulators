#ifndef INITHYDROCARBONSTATE_HPP
#define INITHYDROCARBONSTATE_HPP

#include "opm/core/simulator/BlackoilState.hpp"

namespace Opm
{

void initHydroCarbonState(BlackoilState& state, const PhaseUsage& pu, const int num_cells) {
    enum { Oil = BlackoilPhases::Liquid, Gas = BlackoilPhases::Vapour, Water = BlackoilPhases::Aqua };
    // hydrocarbonstate is only used when gas and oil is present
    assert(pu.phase_used[Oil]);
    std::vector<HydroCarbonState>& hydroCarbonState = state.hydroCarbonState();
    hydroCarbonState.resize(num_cells);
    if (!pu.phase_used[Gas]) {
        // hydroCarbonState should only be used when oil and gas is present. Return OilOnly to avoid potential trouble.
        std::fill(hydroCarbonState.begin(), hydroCarbonState.end(), HydroCarbonState::OilOnly);
        return;
    }
    const int np = pu.num_phases;
    std::fill(hydroCarbonState.begin(), hydroCarbonState.end(), HydroCarbonState::GasAndOil);

    // set hydrocarbon state
    const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
    const std::vector<double>& saturation = state.saturation();
    for (int c = 0; c < num_cells; ++c) {
        if (pu.phase_used[Water]) {
            if ( saturation[c*np + pu.phase_pos[ Water ]] > (1.0 - epsilon)) {
                continue; // cases (almost) filled with water is treated as GasAndOil case;
            }
        }
        if ( saturation[c*np + pu.phase_pos[ Gas ]] == 0.0) {
            hydroCarbonState[c] = HydroCarbonState::OilOnly;
            continue;
        }
        if ( saturation[c*np + pu.phase_pos[ Oil ]] == 0.0) {
            hydroCarbonState[c] = HydroCarbonState::GasOnly;
        }
    }
}


} // namespace Opm
#endif // INITHYDROCARBONSTATE_HPP
