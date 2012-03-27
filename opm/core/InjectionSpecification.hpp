#ifndef OPM_INJECTORSPECIFICATION_HPP
#define	OPM_INJECTORSPECIFICATION_HPP

namespace Opm {

class InjectionSpecification {
public:
    enum Component {
        GAS, OIL, WATER
    };
    
    enum ControlMode {
        NONE, RATE, REIN, RESV, VREP, WGRA, FLD
    };
    
    InjectionSpecification();
    InjectionSpecification(const InjectionSpecification& orig);
    virtual ~InjectionSpecification();
    
    Component component();
    void set_component(Component comp);
    
    ControlMode control_mode();
    void set_control_mode(ControlMode mode);
    
    /// \returns 0 if no limit, else the target/limit of the surface injection 
    ///          rate.
    double surface_injection_target();
    void set_surface_injection_target(double target);
    
    /// \returns 0 if no limit, else the target/limit of the reinjection fraction
    double reinjection_fraction_target();
    void set_reinjection_fraction_target(double target);
    
    /// \returns 0 if no limit, else the target/limit of the BHP
    double BHP_target();
    void set_BHP_target(double target);
private:
    Component component_;
    ControlMode control_mode_;
    double surface_injection_target_;
    double reinjection_fraction_target_;
    double BHP_target_;
};
}
#endif	/* OPM_INJECTORSPECIFICATION_HPP */

