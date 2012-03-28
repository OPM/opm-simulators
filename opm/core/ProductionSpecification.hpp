#ifndef OPM_PRODUCTIONSPECIFICATION_HPP
#define	OPM_PRODUCTIONSPECIFICATION_HPP

#include <opm/core/newwells.h>

namespace Opm
{

    struct ProductionSpecification
    {

        enum ControlMode
        {
            NONE_CM, ORAT, WRAT, REIN, RESV, VREP, WGRA, FLD
        };

        enum Procedure
        {
            NONE_P, RATE, WELL
        };

        ProductionSpecification();
 
        surface_component component_;
        ControlMode control_mode_;
        Procedure procedure_;

        double oil_production_target_;
        double water_production_target_;
        double liquid_production_target_;
        double BHP_target_;


    };
}

#endif	/* OPM_PRODUCTIONSPECIFICATION_HPP */

