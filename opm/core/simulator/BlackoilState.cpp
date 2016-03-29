#include "BlackoilState.hpp"
#include <opm/common/util/numeric/cmp.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>


using namespace Opm;


const std::string BlackoilState::GASOILRATIO = "GASOILRATIO";
const std::string BlackoilState::RV = "RV";
const std::string BlackoilState::SURFACEVOL = "SURFACEVOL";


BlackoilState::BlackoilState( size_t num_cells , size_t num_faces , size_t num_phases)
    : SimulationDataContainer( num_cells , num_faces , num_phases)
{
    registerCellData( GASOILRATIO , 1 );
    registerCellData( RV, 1 );
    registerCellData( SURFACEVOL, num_phases );
    setBlackoilStateReferencePointers();
}

BlackoilState::BlackoilState( const BlackoilState& other )
    : SimulationDataContainer(other)
{
    setBlackoilStateReferencePointers();
}

BlackoilState& BlackoilState::operator=( const BlackoilState& other )
{
    SimulationDataContainer::operator=(other);
    setBlackoilStateReferencePointers();
    return *this;
}

void BlackoilState::setBlackoilStateReferencePointers()
{
    // This sets the reference pointers for the fast
    // accessors, the fields must have been created first.
    gasoilratio_ref_ = &getCellData(GASOILRATIO);
    rv_ref_          = &getCellData(RV);
    surfacevol_ref_  = &getCellData(SURFACEVOL);
}
