#include "BlackoilState.hpp"
#include <opm/core/props/BlackoilPropertiesInterface.hpp>

using namespace Opm;

void
BlackoilState::init(const UnstructuredGrid& g, int num_phases) {
   SimulatorState::init(g, num_phases);
   gor_.resize(g.number_of_cells, 0.) ;
   rv_.resize(g.number_of_cells,0.);
   // surfvol_ intentionally empty, left to initBlackoilSurfvol
}

/// Set the first saturation to either its min or max value in
/// the indicated cells. The second saturation value s2 is set
/// to (1.0 - s1) for each cell. Any further saturation values
/// are unchanged.
void
BlackoilState::setFirstSat(const std::vector<int>& cells,
                           const Opm::BlackoilPropertiesInterface& props,
                           ExtremalSat es) {
    SimulatorState::setFirstSat(cells, props, es);
}

bool
BlackoilState::equals(const SimulatorState& other,
                      double epsilon) const {
    const BlackoilState* that = dynamic_cast <const BlackoilState*> (&other);
    bool equal = that != 0;
    equal = equal && SimulatorState::equals (other, epsilon);
    equal = equal && SimulatorState::vectorApproxEqual(this->surfacevol(),
                                                       that->surfacevol(),
                                                       epsilon);
    equal = equal && SimulatorState::vectorApproxEqual(this->gasoilratio(),
                                                       that->gasoilratio(),
                                                       epsilon);
    return equal;
}
