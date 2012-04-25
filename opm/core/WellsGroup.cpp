/* 
 * File:   WellsGroup.cpp
 * Author: kjetilo
 * 
 * Created on March 27, 2012, 9:27 AM
 */

#include <opm/core/WellsGroup.hpp>
#include <cmath>
#include <opm/core/newwells.h>

namespace Opm
{

    WellPhasesSummed::WellPhasesSummed() 
    : bhp_sum(0.0), rate_sum(0.0)
    {
        
    }
    
    void WellPhasesSummed::operator+=(const WellPhasesSummed& other) 
    {
        rate_sum += other.rate_sum;
        bhp_sum += other.bhp_sum;
    }
    
    WellsGroupInterface::WellsGroupInterface(const std::string& myname,
            ProductionSpecification prod_spec,
            InjectionSpecification inje_spec)
    : parent_(NULL),
    name_(myname),
    production_specification_(prod_spec),
    injection_specification_(inje_spec)
    {
    }

    WellsGroupInterface::~WellsGroupInterface()
    {
    }

    const WellsGroupInterface* WellsGroupInterface::getParent() const
    {
        return parent_;
    }
    const std::string& WellsGroupInterface::name()
    {
        return name_;
    }

    WellsGroup::WellsGroup(const std::string& myname,
            ProductionSpecification prod_spec,
            InjectionSpecification inj_spec)
    : WellsGroupInterface(myname, prod_spec, inj_spec)
    {
    }

    bool WellsGroupInterface::isLeafNode() const
    {
        return false;
    }

    void WellsGroupInterface::setParent(WellsGroupInterface* parent)
    {
        parent_ = parent;
    }

    const ProductionSpecification& WellsGroupInterface::prodSpec() const
    {
        return production_specification_;
    }

    /// Injection specifications for the well or well group.

    const InjectionSpecification& WellsGroupInterface::injSpec() const
    {
        return injection_specification_;
    }

    /// Production specifications for the well or well group.

    ProductionSpecification& WellsGroupInterface::prodSpec()
    {
        return production_specification_;
    }

    /// Injection specifications for the well or well group.

    InjectionSpecification& WellsGroupInterface::injSpec()
    {
        return injection_specification_;
    }

    WellsGroupInterface* WellsGroup::findGroup(const std::string& name_of_node)
    {
        if (name() == name_of_node) {
            return this;
        } else {
            for (size_t i = 0; i < children_.size(); i++) {
                WellsGroupInterface* result = children_[i]->findGroup(name_of_node);
                if (result) {
                    return result;
                }
            }

            // Not found in this node.
            return NULL;
        }
    }

    
    void WellsGroup::calculateGuideRates()
    {
        double guide_rate_sum = 0.0;
        for(size_t i = 0; i < children_.size(); i++) {
            if(children_[i]->isLeafNode()) {
                guide_rate_sum += children_[i]->prodSpec().guide_rate_;
            }
            else
            {
                children_[i]->calculateGuideRates();
            }
        }
        if(guide_rate_sum != 0.0) {
            for(size_t i = 0; i < children_.size(); i++) {
                children_[i]->prodSpec().guide_rate_ /= guide_rate_sum;
            }
        }
    }
    
    void WellsGroup::applyControl(const WellControlType type) 
    {
        for (size_t i = 0; i < children_.size(); ++i) {
            children_[i]->applyControl(type);
        }
    }

    
    bool WellsGroup::conditionsMet(const std::vector<double>& well_bhp, 
                                   const std::vector<double>& well_rate, 
                                   WellPhasesSummed& summed_phases, 
                                   const double epsilon)
    {
        WellPhasesSummed child_phases_summed;
        for(size_t i = 0; i < children_.size(); ++i) {
            WellPhasesSummed current_child_phases_summed;
            if(!children_[i]->conditionsMet(well_bhp, well_rate, 
                                            current_child_phases_summed, epsilon)) {
                return false;
            }
            child_phases_summed += current_child_phases_summed;
        }
        

        double bhp_target = std::min(injSpec().BHP_limit_, prodSpec().BHP_limit_);
        double rate_target = std::min(injSpec().fluid_volume_max_rate_, 
                                      prodSpec().fluid_volume_max_rate_);
        
        double bhp_sum = child_phases_summed.bhp_sum;
        double rate_sum = child_phases_summed.rate_sum;
        if (bhp_sum - bhp_target > epsilon) {
            std::cout << "BHP not met" << std::endl;
            std::cout << "BHP limit was " << bhp_target << std::endl;
            std::cout << "Actual bhp was " << bhp_sum << std::endl;
            
            switch(prodSpec().procedure_) {
            case ProductionSpecification::WELL:
                getWorstOffending(well_bhp).first->shutWell();
                return false;
                break;
            case ProductionSpecification::RATE:
                applyControl(BHP);
                return false;
                break;
            default:
                // Nothing do to;
                break;
            }
        }
        if(rate_sum - rate_target > epsilon) {
            std::cout << "well_rate not met" << std::endl;
            std::cout << "target = " << rate_target 
                      << ", well_rate[index_of_well] = " 
                      << rate_sum << std::endl;
            std::cout << "Group name = " << name() << std::endl;
            
            switch(prodSpec().procedure_) {
            case ProductionSpecification::WELL:
                getWorstOffending(well_rate).first->shutWell();
                return false;
                break;
            case ProductionSpecification::RATE:
                applyControl(RATE);
                return false;
                break;
            default:
                // Nothing do to;
                break;
            }
        }
        
        summed_phases += child_phases_summed;
        return true;
    }

    void WellsGroup::addChild(std::tr1::shared_ptr<WellsGroupInterface> child)
    {
        children_.push_back(child);
    }

    
    int WellsGroup::numberOfLeafNodes() {
        // This could probably use some caching, but seeing as how the number of 
        // wells is relatively small, we'll do without for now.
        int sum = 0;
        
        for(size_t i = 0; i < children_.size(); i++) {
            sum += children_[i]->numberOfLeafNodes();
        }
        
        return sum;
    }
    
    std::pair<WellNode*, double> WellsGroup::getWorstOffending(const std::vector<double>& values) 
    {
        std::pair<WellNode*, double> max;
        for (size_t i = 0; i < children_.size(); i++) {
            std::pair<WellNode*, double> child_max = children_[i]->getWorstOffending(values);
            if (i == 0 || max.second < child_max.second) {
                max = child_max;
            }
        }
        return max;
    }
    
    WellNode::WellNode(const std::string& myname,
            ProductionSpecification prod_spec,
            InjectionSpecification inj_spec)
    : WellsGroupInterface(myname, prod_spec, inj_spec)
    {
    }

    bool WellNode::conditionsMet(const std::vector<double>& well_bhp, 
                                 const std::vector<double>& well_rate, 
                                 WellPhasesSummed& summed_phases, 
                                 const double epsilon)
    {
        
        
        // Check for self:
        if (wells_->type[self_index_] == PRODUCER) {
            double bhp_diff = well_bhp[self_index_] - prodSpec().BHP_limit_;
            double rate_diff = well_rate[self_index_] - prodSpec().fluid_volume_max_rate_;
            
            if (bhp_diff > epsilon) {
                
                std::cout << "BHP exceeded, bhp_diff = " << bhp_diff << std::endl;
                std::cout << "BHP_limit = " << prodSpec().BHP_limit_ << std::endl;
                std::cout << "BHP = " << well_bhp[self_index_] << std::endl; 
                shutWell();
                return false;
            }
            
            if (rate_diff > epsilon) {
                std::cout << "Rate exceeded, rate_diff = " << rate_diff << std::endl;
                shutWell();
                return false;
            }
        } else {
            double bhp_diff = well_bhp[self_index_] - injSpec().BHP_limit_;
            double rate_diff = well_rate[self_index_] - injSpec().fluid_volume_max_rate_;
            
            
           if (bhp_diff > epsilon) {
               std::cout << "BHP exceeded, bhp_diff = " << bhp_diff<<std::endl;
               shutWell();
               return false;
           }
           if (rate_diff > epsilon) {
               std::cout << "Flow diff exceeded, flow_diff = " << rate_diff << std::endl;
               shutWell();
               return false;
           }
            
        }
        
        summed_phases.bhp_sum = well_bhp[self_index_];
        summed_phases.rate_sum = well_rate[self_index_];
        return true;
    }

    WellsGroupInterface* WellNode::findGroup(const std::string& name_of_node)
    {
        if (name() == name_of_node) {
            return this;
        } else {
            return NULL;
        }

    }

    bool WellNode::isLeafNode() const
    {
        return true;
    }

    void WellNode::setWellsPointer(Wells* wells, int self_index)
    {
        wells_ = wells;
        self_index_ = self_index;
    }
    
    void WellNode::calculateGuideRates()
    {
        // Empty
    }
    
    int WellNode::numberOfLeafNodes() 
    {
        return 1;
    }
    
    void WellNode::shutWell() 
    {
        wells_->ctrls[self_index_]->target[0] = 0.0;
    }
    
    std::pair<WellNode*, double> WellNode::getWorstOffending(const std::vector<double>& values) {
        return std::make_pair<WellNode*, double>(this, values[self_index_]);
    }
    
    void WellNode::applyControl(const WellControlType type)
    {
        wells_->ctrls[self_index_]->type[0] = type;
        double target = 0.0;
        switch(type) {
        case BHP:
            if(wells_->type[self_index_] == INJECTOR) {
                target = injSpec().BHP_limit_;
            }
            else {
                target = prodSpec().BHP_limit_;
            }
            break;
        case RATE:
            if(wells_->type[self_index_] == INJECTOR) {
                target = injSpec().fluid_volume_max_rate_;
            }
            else {
                target = prodSpec().fluid_volume_max_rate_;
            }
            break;
        }
        wells_->ctrls[self_index_]->target[0] = target;
    }

    namespace
    {

        SurfaceComponent toSurfaceComponent(std::string type)
        {
            if (type == "OIL") {
                return OIL;
            }
            if (type == "WATER") {
                return WATER;
            }
            if (type == "GAS") {
                return GAS;
            }
            THROW("Unknown type " << type << ", could not convert to SurfaceComponent");
        }

        InjectionSpecification::ControlMode toInjectionControlMode(std::string type)
        {
            if (type == "NONE") {
                return InjectionSpecification::NONE;
            }

            if (type == "ORAT") {
                return InjectionSpecification::ORAT;
            }
            if (type == "REIN") {
                return InjectionSpecification::REIN;
            }
            if (type == "RESV") {
                return InjectionSpecification::RESV;
            }
            if (type == "VREP") {
                return InjectionSpecification::VREP;
            }
            if (type == "WGRA") {
                return InjectionSpecification::WGRA;
            }
            if (type == "FLD") {
                return InjectionSpecification::FLD;
            }
            if (type == "GRUP") {
                return InjectionSpecification::GRUP;
            }


            THROW("Unknown type " << type << ", could not convert to ControlMode.");
        }

        ProductionSpecification::ControlMode toProductionControlMode(std::string type)
        {
            if (type == "NONE") {
                return ProductionSpecification::NONE_CM;
            }
            if (type == "ORAT") {
                return ProductionSpecification::ORAT;

            }
            if (type == "LRAT") {
                return ProductionSpecification::LRAT;
            }
            if (type == "REIN") {
                return ProductionSpecification::REIN;
            }
            if (type == "RESV") {
                return ProductionSpecification::RESV;
            }
            if (type == "VREP") {
                return ProductionSpecification::VREP;
            }
            if (type == "WGRA") {
                return ProductionSpecification::WGRA;
            }
            if (type == "FLD") {
                return ProductionSpecification::FLD;
            }
            if (type == "GRUP") {
                return ProductionSpecification::GRUP;
            }
            if (type == "BHP") {
                return ProductionSpecification::BHP;
            }

            THROW("Unknown type " << type << ", could not convert to ControlMode.");
        }

        ProductionSpecification::Procedure toProductionProcedure(std::string type)
        {
            if (type == "NONE") {
                return ProductionSpecification::NONE_P;
            }
            if (type == "RATE") {
                return ProductionSpecification::RATE;
            }
            if (type == "WELL") {
                return ProductionSpecification::WELL;
            }


            THROW("Unknown type " << type << ", could not convert to ControlMode.");
        }
    } // anonymous namespace

    std::tr1::shared_ptr<WellsGroupInterface> createWellsGroup(const std::string& name, const EclipseGridParser& deck)
    {

        std::tr1::shared_ptr<WellsGroupInterface> return_value;
        // First we need to determine whether it's a group or just a well:
        bool isWell = false;
        if (deck.hasField("WELSPECS")) {
            WELSPECS wspecs = deck.getWELSPECS();
            for (size_t i = 0; i < wspecs.welspecs.size(); i++) {
                if (wspecs.welspecs[i].name_ == name) {
                    isWell = true;
                    break;
                }
            }
        }
        // For now, assume that if it isn't a well, it's a group

        if (isWell) {

            InjectionSpecification injection_specification;
            if (deck.hasField("WCONINJE")) {
                WCONINJE wconinje = deck.getWCONINJE();
                for (size_t i = 0; i < wconinje.wconinje.size(); i++) {
                    if (wconinje.wconinje[i].well_ == name) {
                        WconinjeLine line = wconinje.wconinje[i];
                        injection_specification.BHP_limit_ = line.BHP_limit_;
                        injection_specification.injector_type_ = toSurfaceComponent(line.injector_type_);
                        injection_specification.control_mode_ = toInjectionControlMode(line.control_mode_);
                        injection_specification.surface_flow_max_rate_ = line.surface_flow_max_rate_;
                        injection_specification.fluid_volume_max_rate_ = line.fluid_volume_max_rate_;
                    }
                }
            }

            ProductionSpecification production_specification;
            if (deck.hasField("WCONPROD")) {
                WCONPROD wconprod = deck.getWCONPROD();
                std::cout << wconprod.wconprod.size() << std::endl;
                for (size_t i = 0; i < wconprod.wconprod.size(); i++) {
                    if (wconprod.wconprod[i].well_ == name) {
                        WconprodLine line = wconprod.wconprod[i];
                        production_specification.BHP_limit_ = line.BHP_limit_;
                        production_specification.fluid_volume_max_rate_ = line.fluid_volume_max_rate_;
                        production_specification.oil_max_rate_ = line.oil_max_rate_;
                        production_specification.control_mode_ = toProductionControlMode(line.control_mode_);
                        production_specification.water_production_target_ = line.water_max_rate_;
                    }
                }
            }
            return_value.reset(new WellNode(name, production_specification, injection_specification));
        } else {
            InjectionSpecification injection_specification;
            if (deck.hasField("GCONINJE")) {
                GCONINJE gconinje = deck.getGCONINJE();
                for (size_t i = 0; i < gconinje.gconinje.size(); i++) {
                    if (gconinje.gconinje[i].group_ == name) {
                        GconinjeLine line = gconinje.gconinje[i];
                        injection_specification.injector_type_ = toSurfaceComponent(line.injector_type_);
                        injection_specification.control_mode_ = toInjectionControlMode(line.control_mode_);
                        injection_specification.surface_flow_max_rate_ = line.surface_flow_max_rate_;
                        injection_specification.fluid_volume_max_rate_ = line.resv_flow_max_rate_;
                    }
                }
            }

            ProductionSpecification production_specification;
            if (deck.hasField("GCONPROD")) {
                std::cout << "Searching in gconprod " << std::endl;
                std::cout << "name= " << name << std::endl;
                GCONPROD gconprod = deck.getGCONPROD();
                for (size_t i = 0; i < gconprod.gconprod.size(); i++) {
                    if (gconprod.gconprod[i].group_ == name) {
                        GconprodLine line = gconprod.gconprod[i];
                        production_specification.oil_max_rate_ = line.oil_max_rate_;
                        std::cout << "control_mode = " << line.control_mode_ << std::endl;
                        production_specification.control_mode_ = toProductionControlMode(line.control_mode_);
                        production_specification.water_production_target_ = line.water_max_rate_;
                        production_specification.gas_max_rate_ = line.gas_max_rate_;
                        production_specification.liquid_max_rate_ = line.liquid_max_rate_;
                        production_specification.procedure_ = toProductionProcedure(line.procedure_);
                    }
                }
            }

            return_value.reset(new WellsGroup(name, production_specification, injection_specification));
        }

        return return_value;
    }
}
