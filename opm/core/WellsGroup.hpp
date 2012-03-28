#ifndef OPM_WELLSGROUP_HPP
#define	OPM_WELLSGROUP_HPP

#include <opm/core/InjectionSpecification.hpp>
#include <opm/core/ProductionSpecification.hpp>
#include <string>

namespace Opm
{

    class WellsGroupInterface
    {
    public:
        WellsGroupInterface(const std::string& name);
        virtual ~WellsGroupInterface();

        const std::string& name();
        const ProductionSpecification& prodSpec() const;
        const InjectionSpecification& injSpec() const;

    private:
        std::string name_;
        ProductionSpecification production_specification_;
        InjectionSpecification injection_specification_;

    };

    class WellsGroup : public WellsGroupInterface
    {
    };

}
#endif	/* OPM_WELLSGROUP_HPP */

