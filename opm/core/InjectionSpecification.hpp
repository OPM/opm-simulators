#ifndef OPM_INJECTORSPECIFICATION_HPP
#define	OPM_INJECTORSPECIFICATION_HPP

#include <opm/core/newwells.h>

namespace Opm
{

    struct InjectionSpecification
    {

        enum ControlMode
        {
            NONE, ORAT, REIN, RESV, VREP, WGRA, FLD, GRUP
        };

        enum InjectorType
        {
            WATER, OIL, GAS
        };

        InjectionSpecification();

        InjectorType injector_type_;
        ControlMode control_mode_;
        double surface_flow_max_rate_;
        double reinjection_fraction_target_;
        double fluid_volume_max_rate_;
        double BHP_limit_;
    };

}

#endif	/* OPM_INJECTORSPECIFICATION_HPP */

