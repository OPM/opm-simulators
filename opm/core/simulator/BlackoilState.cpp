#include "BlackoilState.hpp"
#include <opm/core/props/BlackoilPropertiesInterface.hpp>

using namespace Opm;

void
BlackoilState::init(int number_of_cells, int number_of_phases, int num_phases) {
    SimulatorState::init(number_of_cells, number_of_phases, num_phases);
   gor_.resize(number_of_cells, 0.) ;
   rv_.resize(number_of_cells,0.);
   // surfvol_ intentionally empty, left to initBlackoilSurfvol
}

void
BlackoilState::init(const UnstructuredGrid& g, int num_phases)
{
    init(g.number_of_cells, g.number_of_faces, num_phases);
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
