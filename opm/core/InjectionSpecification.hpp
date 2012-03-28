#ifndef OPM_INJECTORSPECIFICATION_HPP
#define	OPM_INJECTORSPECIFICATION_HPP

#include <opm/core/newwells.h>
namespace Opm
{

    struct InjectionSpecification
    {

        enum ControlMode
        {
            NONE, RATE, REIN, RESV, VREP, WGRA, FLD
        };

        InjectionSpecification();

        surface_component component_;
        ControlMode control_mode_;
        double surface_injection_target_;
        double reinjection_fraction_target_;
        double BHP_target_;
    };
}
#endif	/* OPM_INJECTORSPECIFICATION_HPP */

