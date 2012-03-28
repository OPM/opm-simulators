#ifndef OPM_PRODUCTIONSPECIFICATION_HPP
#define	OPM_PRODUCTIONSPECIFICATION_HPP

namespace Opm {
class ProductionSpecification {
public:
    enum Component {
        GAS, OIL, WATER
    };
    
    enum ControlMode {
        NONE_CM, ORAT, WRAT, REIN, RESV, VREP, WGRA, FLD
    };
    
    enum Procedure {
        WELL, RATE, NONE_P
    };
    
    ProductionSpecification();
    virtual ~ProductionSpecification();
    
    
    void set_BHP_target(double BHP_target_);
    double get_BHP_target() const;
    void set_liquid_production_target(double liquid_production_target_);
    double get_liquid_production_target() const;
    void set_water_production_target(double water_production_target_);
    double get_water_production_target() const;
    void set_oil_production_target(double oil_production_target_);
    double get_oil_production_target() const;
    void set_procedure(Procedure procedure_);
    Procedure get_procedure() const;
    void set_control_mode(ControlMode control_mode_);
    ControlMode get_control_mode() const;
    void set_component(Component component_);
    Component get_component() const;
private:
    Component component_;
    ControlMode control_mode_;
    Procedure procedure_;

    double oil_production_target_;
    double water_production_target_;
    double liquid_production_target_;
    double BHP_target_;
    
    
};
}

#endif	/* OPM_PRODUCTIONSPECIFICATION_HPP */

