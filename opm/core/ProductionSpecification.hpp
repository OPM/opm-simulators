#ifndef OPM_PRODUCTIONSPECIFICATION_HPP
#define	OPM_PRODUCTIONSPECIFICATION_HPP

#include <opm/core/newwells.h>

namespace Opm
{

    struct ProductionSpecification
    {

        enum ControlMode
        {
            NONE_CM, ORAT, WRAT, LRAT, REIN, RESV, VREP, WGRA, FLD, GRUP
        };

        enum Procedure
        {
            NONE_P, RATE, WELL
        };
        
        enum GuideRateType
        {
            OIL, RAT, NONE_GRT
        };

        ProductionSpecification();
 
        ControlMode control_mode_;
        Procedure procedure_;

        double oil_max_rate_;
        double water_production_target_;
        double fluid_volume_max_rate_;
        double BHP_limit_;
        double gas_max_rate_;
        double liquid_max_rate_;
        double guide_rate_;
        GuideRateType guide_rate_type_;

    };
}

#endif	/* OPM_PRODUCTIONSPECIFICATION_HPP */

