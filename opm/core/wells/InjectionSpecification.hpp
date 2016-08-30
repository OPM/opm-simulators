#ifndef OPM_INJECTORSPECIFICATION_HPP
#define	OPM_INJECTORSPECIFICATION_HPP

#include <opm/core/wells.h>
#include <string>

namespace Opm
{

    struct InjectionSpecification
    {

        enum ControlMode
        {
            NONE, RATE, RESV, BHP, THP, REIN, VREP, GRUP, FLD
        };

        enum InjectorType
        {
            WATER, OIL, GAS
        };

        enum GuideRateType
        {
            RAT, NONE_GRT
        };

        InjectionSpecification();
        static std::string toString(const ControlMode& mode);
        static std::string toString(const InjectorType& type);
        static std::string toString(const GuideRateType& type);
        InjectorType injector_type_;
        ControlMode control_mode_;
        double surface_flow_max_rate_;
        double reservoir_flow_max_rate_;
        double BHP_limit_;
        double reinjection_fraction_target_;
        double voidage_replacment_fraction_;
        double guide_rate_;
        GuideRateType guide_rate_type_;
    };

}

#endif	/* OPM_INJECTORSPECIFICATION_HPP */

