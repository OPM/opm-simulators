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
}
