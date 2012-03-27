#ifndef OPM_WELLSGROUP_HPP
#define	OPM_WELLSGROUP_HPP
#include <opm/core/InjectionSpecification.hpp>
#include <opm/core/ProductionSpecification.hpp>


namespace Opm {
class AbstractWellsGroup {
public:
    AbstractWellsGroup(const std::string& name);
    virtual ~AbstractWellsGroup();
    
    const std::string& get_name();
    const ProductionSpecification& get_production_specification() const;
    const InjectionSpecification& get_injection_specification() const;
    
private:
    std::string name_;
    ProductionSpecification production_specification_;
    InjectionSpecification injection_specification_;

};

class WellsGroup : public AbstractWellsGroup {
    
};

}
#endif	/* OPM_WELLSGROUP_HPP */

