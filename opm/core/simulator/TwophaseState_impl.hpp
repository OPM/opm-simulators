#ifndef OPM_TWOPHASESTATE_HEADER_INCLUDED
#error Do not include this file directly!
#endif

namespace Opm {

inline void
TwophaseState::setFirstSat(const std::vector<int>& cells,
                           const Opm::IncompPropertiesInterface& props,
                           ExtremalSat es) {
    SimulatorState::setFirstSat(cells, props, es);
}

inline bool
TwophaseState::equals (const SimulatorState& other,
                       double epsilon) const {
    return dynamic_cast <const TwophaseState*>(&other)
            ? SimulatorState::equals (other, epsilon)
            : false;
}

} // namespace Opm
