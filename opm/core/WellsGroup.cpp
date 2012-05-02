/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/core/WellsGroup.hpp>
#include <cmath>
#include <opm/core/newwells.h>
#include <opm/core/fluid/blackoil/phaseUsageFromDeck.hpp>

namespace Opm
{

    // ==========   WellPhasesSummed methods   ===========

    WellPhasesSummed::WellPhasesSummed() 
    {
        for (int i = 0; i < 3; ++i) {
            res_inj_rates[i] = 0.0;
            res_prod_rates[i] = 0.0;
            surf_inj_rates[i] = 0.0;
            surf_prod_rates[i] = 0.0;
        }
    }

    void WellPhasesSummed::operator+=(const WellPhasesSummed& other) 
    {
        for (int i = 0; i < 3; ++i) {
            res_inj_rates[i] += other.res_inj_rates[i];
            res_prod_rates[i] += other.res_prod_rates[i];
            surf_inj_rates[i] += other.surf_inj_rates[i];
            surf_prod_rates[i] += other.surf_prod_rates[i];
        }
    }

    // ==========   WellsGroupInterface methods   ===========


    WellsGroupInterface::WellsGroupInterface(const std::string& myname,
                                             const ProductionSpecification& prod_spec,
                                             const InjectionSpecification& inje_spec,
                                             const PhaseUsage& phase_usage)
        : parent_(NULL),
          name_(myname),
          production_specification_(prod_spec),
          injection_specification_(inje_spec),
          phase_usage_(phase_usage)
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

    /// Calculates the correct rate for the given ProductionSpecification::ControlMode
    double WellsGroupInterface::rateByMode(const double* res_rates, 
                                           const double* surf_rates,
                                           const ProductionSpecification::ControlMode mode)
    {
        switch (mode) {
        case ProductionSpecification::ORAT:
            return surf_rates[phaseUsage().phase_pos[BlackoilPhases::Liquid]];
        case ProductionSpecification::WRAT:
            return surf_rates[phaseUsage().phase_pos[BlackoilPhases::Aqua]];
        case ProductionSpecification::GRAT:
            return surf_rates[phaseUsage().phase_pos[BlackoilPhases::Vapour]];
        case ProductionSpecification::LRAT:
            return surf_rates[phaseUsage().phase_pos[BlackoilPhases::Liquid]] 
                + surf_rates[phaseUsage().phase_pos[BlackoilPhases::Aqua]];
        case ProductionSpecification::RESV:
            {
                double tot_rate = 0.0;
                for (int phase = 0; phase < phaseUsage().num_phases; ++phase) {
                    tot_rate += res_rates[phase];
                }
                return tot_rate;
            }
        default:
            THROW("No rate associated with production control mode" << mode);
        }
    }

    /// Calculates the correct rate for the given InjectionSpecification::ControlMode
    double WellsGroupInterface::rateByMode(const double* res_rates, 
                                           const double* surf_rates,
                                           const InjectionSpecification::ControlMode mode)
    {
        const double* rates = 0;
        switch (mode) {
        case InjectionSpecification::RATE:
            rates = surf_rates;
            break;
        case InjectionSpecification::RESV:
            rates = res_rates;
            break;
        default:
            THROW("No rate associated with injection control mode" << mode);
        }
        double tot_rate = 0.0;
        for (int phase = 0; phase < phaseUsage().num_phases; ++phase) {
            tot_rate += rates[phase];
        }
        return tot_rate;
    }


   

    // ==============   WellsGroup members =============

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


    WellsGroup::WellsGroup(const std::string& myname,
                           const ProductionSpecification& prod_spec,
                           const InjectionSpecification& inj_spec,
                           const PhaseUsage& phase_usage)
        : WellsGroupInterface(myname, prod_spec, inj_spec, phase_usage)
    {
    }


    void WellsGroup::calculateGuideRates()
    {
        double inj_guide_rate_sum = 0.0;
        double prod_guide_rate_sum = 0.0;
        for (size_t i = 0; i < children_.size(); i++) {
            children_[i]->calculateGuideRates();
            inj_guide_rate_sum += children_[i]->injSpec().guide_rate_;
            prod_guide_rate_sum += children_[i]->prodSpec().guide_rate_;
        }
        injSpec().guide_rate_ = inj_guide_rate_sum;
        prodSpec().guide_rate_ = prod_guide_rate_sum;
    }
    


    /// Sets the current active control to the provided one for all injectors within the group.
    /// After this call, the combined rate (which rate depending on control_mode) of the group
    /// shall be equal to target.
    void WellsGroup::applyInjGroupControl(const InjectionSpecification::ControlMode control_mode,
                                          const double target)
    {
        for (size_t i = 0; i < children_.size(); ++i) {
            const double child_target = target * children_[i]->injSpec().guide_rate_/injSpec().guide_rate_;
            children_[i]->applyInjGroupControl(control_mode, child_target);
        }
        injSpec().control_mode_ = InjectionSpecification::FLD;
    }

    /// Sets the current active control to the provided one for all producers within the group.
    /// After this call, the combined rate (which rate depending on control_mode) of the group
    /// shall be equal to target.
    void WellsGroup::applyProdGroupControl(const ProductionSpecification::ControlMode control_mode,
                                           const double target)
    {
        for (size_t i = 0; i < children_.size(); ++i) {
            const double child_target = target * children_[i]->prodSpec().guide_rate_/prodSpec().guide_rate_;
            children_[i]->applyProdGroupControl(control_mode, child_target);
        }
        prodSpec().control_mode_ = ProductionSpecification::FLD;
    }


    bool WellsGroup::conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_reservoirrates_phase,
                                   const std::vector<double>& well_surfacerates_phase,
                                   WellPhasesSummed& summed_phases)
    {
        // Check children's constraints recursively.
        WellPhasesSummed child_phases_summed;
        for (size_t i = 0; i < children_.size(); ++i) {
            WellPhasesSummed current_child_phases_summed;
            if (!children_[i]->conditionsMet(well_bhp,
                                             well_reservoirrates_phase,
                                             well_surfacerates_phase,
                                             current_child_phases_summed)) {
                return false;
            }
            child_phases_summed += current_child_phases_summed;
        }

        const int np = phaseUsage().num_phases;

        // Injection constraints.
        // RATE
        if (injSpec().control_mode_ != InjectionSpecification::RATE) {
            const double target_rate = injSpec().surface_flow_max_rate_;
            if (target_rate >= 0.0) {
                double my_rate = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    my_rate += child_phases_summed.surf_inj_rates[phase];
                }
                if (my_rate > target_rate) {
                    std::cout << "Group RATE target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyInjGroupControl(InjectionSpecification::RATE, target_rate);
                    injSpec().control_mode_ = InjectionSpecification::RATE;
                    return false;
                }
            }
        }
        // RESV
        if (injSpec().control_mode_ != InjectionSpecification::RESV) {
            const double target_rate = injSpec().reservoir_flow_max_rate_;
            if (target_rate >= 0.0) {
                double my_rate = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    my_rate += child_phases_summed.res_inj_rates[phase];
                }
                if (my_rate > target_rate) {
                    std::cout << "Group RESV target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyInjGroupControl(InjectionSpecification::RESV, target_rate);
                    injSpec().control_mode_ = InjectionSpecification::RESV;
                    return false;
                }
            }
        }
        // REIN
        // \TODO: Add support for REIN controls.

        // Production constraints.
        bool prod_restrictions_violated = false;
        ProductionSpecification::ControlMode violated_prod_mode = ProductionSpecification::NONE;

        rateByMode(child_phases_summed.res_prod_rates,
                   child_phases_summed.surf_prod_rates,
                   mode);
        
        // ORAT
        if (prodSpec().control_mode_ != ProductionSpecification::ORAT) {
            const double target_rate = prodSpec().oil_max_rate_;
            if (target_rate >= 0.0) {
                const double my_rate
                    = child_phases_summed.surf_prod_rates[phaseUsage().phase_pos[BlackoilPhases::Liquid]];
                if (std::fabs(my_rate) > target_rate) {
                    std::cout << "Group ORAT target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyProdGroupControl(ProductionSpecification::ORAT, target_rate);
                    prodSpec().control_mode_ = ProductionSpecification::ORAT;
                    return false;
                }
            }
        }

        // WRAT
        if (prodSpec().control_mode_ != ProductionSpecification::WRAT) {
            const double target_rate = prodSpec().water_max_rate_;
            if (target_rate >= 0.0) {
                const double my_rate
                    = child_phases_summed.surf_prod_rates[phaseUsage().phase_pos[BlackoilPhases::Aqua]];
                if (std::fabs(my_rate) > target_rate) {
                    std::cout << "Group WRAT target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyProdGroupControl(ProductionSpecification::WRAT, target_rate);
                    prodSpec().control_mode_ = ProductionSpecification::WRAT;
                    return false;
                }
            }
        }

        // GRAT
        if (prodSpec().control_mode_ != ProductionSpecification::GRAT) {
            const double target_rate = prodSpec().gas_max_rate_;
            if (target_rate >= 0.0) {
                const double my_rate
                    = child_phases_summed.surf_prod_rates[phaseUsage().phase_pos[BlackoilPhases::Vapour]];
                if (std::fabs(my_rate) > target_rate) {
                    std::cout << "Group GRAT target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyProdGroupControl(ProductionSpecification::GRAT, target_rate);
                    prodSpec().control_mode_ = ProductionSpecification::GRAT;
                    return false;
                }
            }
        }

        // LRAT
        if (prodSpec().control_mode_ != ProductionSpecification::LRAT) {
            const double target_rate = prodSpec().liquid_max_rate_;
            if (target_rate >= 0.0) {
                const double my_rate = 
                    = child_phases_summed.surf_prod_rates[phaseUsage().phase_pos[BlackoilPhases::Aqua]]
                    + child_phases_summed.surf_prod_rates[phaseUsage().phase_pos[BlackoilPhases::Liquid]];
                if (std::fabs(my_rate) > target_rate) {
                    std::cout << "Group LRAT target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyProdGroupControl(ProductionSpecification::LRAT, target_rate);
                    prodSpec().control_mode_ = ProductionSpecification::LRAT;
                    return false;
                }
            }
        }

        // RESV
        if (prodSpec().control_mode_ != ProductionSpecification::RESV) {
            const double target_rate = prodSpec().reservoir_flow_max_rate_;
            if (target_rate >= 0.0) {
                double my_rate = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    my_rate += child_phases_summed.res_prod_rates[phase];
                }
                if (std::fabs(my_rate) > target_rate) {
                    std::cout << "Group RESV target not met for group " << name() << std::endl;
                    std::cout << "target = " << target_rate << '\n'
                              << "rate   = " << my_rate << std::endl;
                    applyProdGroupControl(ProductionSpecification::RESV, target_rate);
                    prodSpec().control_mode_ = ProductionSpecification::RESV;
                    return false;
                }
            }
        }


        double rate_target = std::min(std::abs(injSpec().fluid_volume_max_rate_), 
                                      prodSpec().fluid_volume_max_rate_);
        
        double rate_sum = child_phases_summed.rate_sum;
        if (std::abs(rate_sum) - std::abs(rate_target) > epsilon) {
            std::cout << "well_rate not met" << std::endl;
            std::cout << "target = " << rate_target 
                      << ", well_rate[index_of_well] = " 
                      << rate_sum << std::endl;
            std::cout << "Group name = " << name() << std::endl;
            
            switch (prodSpec().procedure_) {
            case ProductionSpecification::WELL:
                getWorstOffending(well_rate).first->shutWell();
                return false;
                break;
            case ProductionSpecification::RATE:
                applyControl(SURFACE_RATE);
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



    // ==============    WellNode members   ============


    
    WellNode::WellNode(const std::string& myname,
                       const ProductionSpecification& prod_spec,
                       const InjectionSpecification& inj_spec,
                       const PhaseUsage& phase_usage)
        : WellsGroupInterface(myname, prod_spec, inj_spec, phase_usage),
          wells_(0),
          self_index_(-1),
          group_control_index_(-1)
    {
    }

    bool WellNode::conditionsMet(const std::vector<double>& well_bhp,
                                 const std::vector<double>& well_reservoirrates_phase,
                                 const std::vector<double>& well_surfacerates_phase,
                                 WellPhasesSummed& summed_phases)
    {
        // Report on our rates.
        const int np = phaseUsage().num_phases;
        for (int phase = 0; phase < np; ++phase) {
            if (wells_->type[self_index_] == INJECTOR) {
                summed_phases.res_inj_rates[phase] = well_reservoirrates_phase[np*self_index_ + phase];
                summed_phases.surf_inj_rates[phase] = well_surfacerates_phase[np*self_index_ + phase];
            } else {
                summed_phases.res_prod_rates[phase] = well_reservoirrates_phase[np*self_index_ + phase];
                summed_phases.surf_prod_rates[phase] = well_surfacerates_phase[np*self_index_ + phase];
            }
        }

        // Check constraints.
        bool is_producer = (wells_->type[self_index_] == PRODUCER);
        const WellControls& ctrls = *wells_->ctrls[self_index_];
        for (int ctrl_index = 0; ctrl_index < ctrls.num; ++ctrl_index) {
            if (ctrl_index == ctrls.current || ctrl_index == group_control_index_) {
                // We do not check constraints that either were used
                // as the active control, or that come from group control.
                continue;
            }
            bool ctrl_violated = false;
            switch (ctrls.type[ctrl_index]) {
            case BHP: {
                const double my_well_bhp = well_bhp[self_index_];
                const double my_target_bhp = ctrls.target[ctrl_index];
                ctrl_violated = is_producer ? (my_target_bhp > my_well_bhp)
                    : (my_target_bhp < my_well_bhp);
                if (ctrl_violated) {
                    std::cout << "BHP limit violated for well " << name() << ":\n";
                    std::cout << "BHP limit = " << my_target_bhp << std::endl;
                    std::cout << "BHP       = " << my_well_bhp << std::endl;
                }
                break;
            }
            case RESERVOIR_RATE: {
                double my_rate = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    my_rate += ctrls.distr[np*ctrl_index + phase]*well_reservoirrates_phase[np*self_index_ + phase];
                }
                const double my_rate_target = ctrls.target[ctrl_index];
                ctrl_violated = std::fabs(my_rate) > std::fabs(my_rate_target);
                if (ctrl_violated) {
                    std::cout << "RESERVOIR_RATE limit violated for well " << name() << ":\n";
                    std::cout << "rate limit = " << my_rate_target;
                    std::cout << "rate       = " << my_rate;
                }
                break;
            }
            case SURFACE_RATE: {
                double my_rate = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    my_rate += ctrls.distr[np*ctrl_index + phase]*well_surfacerates_phase[np*self_index_ + phase];
                }
                const double my_rate_target = ctrls.target[ctrl_index];
                ctrl_violated = std::fabs(my_rate) > std::fabs(my_rate_target);
                if (ctrl_violated) {
                    std::cout << "SURFACE_RATE limit violated for well " << name() << ":\n";
                    std::cout << "rate limit = " << my_rate_target;
                    std::cout << "rate       = " << my_rate;
                }
                break;
            }
            } // end of switch()
            if (ctrl_violated) {
                set_current_control(self_index_, ctrl_index, wells_);
                return false;
            }
        }
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
        bool already_has_group_control = 
            ((wells_->type[self_index_] == INJECTOR) && (injSpec().control_mode_ == InjectionSpecification::GRUP))
            || ((wells_->type[self_index_] == PRODUCER) && (prodSpec().control_mode_ == ProductionSpecification::GRUP));
        if (already_has_group_control) {
            group_control_index_ = wells_->ctrls[self_index_]->num - 1;
        }
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
    
    void WellNode::applyInjGroupControl(const InjectionSpecification::ControlMode control_mode,
                                        const double target)
    {
        if (!wells_->type[self_index_] == INJECTOR) {
            ASSERT(target == 0.0);
            return;
        }

        const double distr[3] = { 1.0, 1.0, 1.0 };
        WellControlType wct;
        switch (control_mode) {
        case InjectionSpecification::RATE:
            wct = SURFACE_RATE;
            break;
        case InjectionSpecification::RESV:
            wct = RESERVOIR_RATE;
            break;
        default:
            THROW("Group injection control mode not handled: " << control_mode);
        }

        if (group_control_index_ < 0) {
            // The well only had its own controls, no group controls.
            append_well_controls(wct, target, distr, self_index_, wells_);
            group_control_index_ = wells_->ctrls[self_index_]->num - 1;
        } else {
            // We will now modify the last control, that
            // "belongs to" the group control.
            const int np = wells_->number_of_phases;
            wells_->ctrls[self_index_]->type[group_control_index_] = wct;
            wells_->ctrls[self_index_]->target[group_control_index_] = target;
            std::copy(distr, distr + np, wells_->ctrls[self_index_]->distr + np*group_control_index_);
        }
        set_current_control(self_index_, group_control_index_, wells_);
    }


    void WellNode::applyProdGroupControl(const ProductionSpecification::ControlMode control_mode,
                                         const double target)
    {
        if (!wells_->type[self_index_] == PRODUCER) {
            ASSERT(target == 0.0);
            return;
        }

        double distr[3] = { 0.0, 0.0, 0.0 };
        const int* phase_pos = phaseUsage().phase_pos;
        const int* phase_used = phaseUsage().phase_used;
        WellControlType wct;
        switch (control_mode) {
        case ProductionSpecification::ORAT:
            wct = SURFACE_RATE;
            if (!phase_used[BlackoilPhases::Liquid]) {
                THROW("Oil phase not active and ORAT control specified.");
            }
            distr[phase_pos[BlackoilPhases::Liquid]] = 1.0;
            break;
        case ProductionSpecification::WRAT:
            wct = SURFACE_RATE;
            if (!phase_used[BlackoilPhases::Aqua]) {
                THROW("Water phase not active and WRAT control specified.");
            }
            distr[phase_pos[BlackoilPhases::Aqua]] = 1.0;
            break;
        case ProductionSpecification::GRAT:
            wct = SURFACE_RATE;
            if (!phase_used[BlackoilPhases::Vapour]) {
                THROW("Gas phase not active and GRAT control specified.");
            }
            distr[phase_pos[BlackoilPhases::Vapour]] = 1.0;
            break;
        case ProductionSpecification::LRAT:
            wct = SURFACE_RATE;
            if (!phase_used[BlackoilPhases::Liquid]) {
                THROW("Oil phase not active and LRAT control specified.");
            }
            if (!phase_used[BlackoilPhases::Aqua]) {
                THROW("Water phase not active and LRAT control specified.");
            }
            distr[phase_pos[BlackoilPhases::Liquid]] = 1.0;
            distr[phase_pos[BlackoilPhases::Aqua]] = 1.0;
            break;
        case ProductionSpecification::RESV:
            distr[0] = distr[1] = distr[2] = 1.0;
            wct = RESERVOIR_RATE;
            break;
        default:
            THROW("Group injection control mode not handled: " << control_mode);
        }

        if (group_control_index_ < 0) {
            // The well only had its own controls, no group controls.
            append_well_controls(wct, target, distr, self_index_, wells_);
            group_control_index_ = wells_->ctrls[self_index_]->num - 1;
        } else {
            // We will now modify the last control, that
            // "belongs to" the group control.
            const int np = wells_->number_of_phases;
            wells_->ctrls[self_index_]->type[group_control_index_] = wct;
            wells_->ctrls[self_index_]->target[group_control_index_] = target;
            std::copy(distr, distr + np, wells_->ctrls[self_index_]->distr + np*group_control_index_);
        }
        set_current_control(self_index_, group_control_index_, wells_);
    }

    namespace
    {

        InjectionSpecification::InjectorType toInjectorType(std::string type)
        {
            if (type == "OIL") {
                return InjectionSpecification::OIL;
            }
            if (type == "WATER") {
                return InjectionSpecification::WATER;
            }
            if (type == "GAS") {
                return InjectionSpecification::GAS;
            }
            THROW("Unknown type " << type << ", could not convert to SurfaceComponent");
        }


#define HANDLE_ICM(x)                           \
        if (type == #x) {                       \
            return InjectionSpecification::x;   \
        }

        InjectionSpecification::ControlMode toInjectionControlMode(std::string type)
        {
            HANDLE_ICM(NONE);
            HANDLE_ICM(RATE);
            HANDLE_ICM(RESV);
            HANDLE_ICM(BHP);
            HANDLE_ICM(THP);
            HANDLE_ICM(REIN);
            HANDLE_ICM(VREP);
            HANDLE_ICM(GRUP);
            HANDLE_ICM(FLD);
            THROW("Unknown type " << type << ", could not convert to InjectionSpecification::ControlMode.");
        }
#undef HANDLE_ICM

#define HANDLE_PCM(x)                           \
        if (type == #x) {                       \
            return ProductionSpecification::x;  \
        }

        ProductionSpecification::ControlMode toProductionControlMode(std::string type)
        {
            HANDLE_PCM(NONE);
            HANDLE_PCM(ORAT);
            HANDLE_PCM(WRAT);
            HANDLE_PCM(GRAT);
            HANDLE_PCM(LRAT);
            HANDLE_PCM(CRAT);
            HANDLE_PCM(RESV);
            HANDLE_PCM(PRBL);
            HANDLE_PCM(BHP);
            HANDLE_PCM(THP);
            HANDLE_PCM(GRUP);
            HANDLE_PCM(FLD);
            THROW("Unknown type " << type << ", could not convert to ProductionSpecification::ControlMode.");
        }
#undef HANDLE_PCM

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

    std::tr1::shared_ptr<WellsGroupInterface> createWellsGroup(const std::string& name,
                                                               const EclipseGridParser& deck)
    {
        PhaseUsage phase_usage = phaseUsageFromDeck(deck);

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
                        injection_specification.injector_type_ = toInjectorType(line.injector_type_);
                        injection_specification.control_mode_ = toInjectionControlMode(line.control_mode_);
                        injection_specification.surface_flow_max_rate_ = line.surface_flow_max_rate_;
                        injection_specification.reservoir_flow_max_rate_ = line.reservoir_flow_max_rate_;
                    }
                }
            } else {
                injection_specification.guide_rate_ = 0.0;
            }

            ProductionSpecification production_specification;
            if (deck.hasField("WCONPROD")) {
                WCONPROD wconprod = deck.getWCONPROD();
                std::cout << wconprod.wconprod.size() << std::endl;
                for (size_t i = 0; i < wconprod.wconprod.size(); i++) {
                    if (wconprod.wconprod[i].well_ == name) {
                        WconprodLine line = wconprod.wconprod[i];
                        production_specification.BHP_limit_ = line.BHP_limit_;
                        production_specification.reservoir_flow_max_rate_ = line.reservoir_flow_max_rate_;
                        production_specification.oil_max_rate_ = line.oil_max_rate_;
                        production_specification.control_mode_ = toProductionControlMode(line.control_mode_);
                        production_specification.water_max_rate_ = line.water_max_rate_;
                    }
                }
            } else {
                production_specification.guide_rate_ = 0.0;
            }
            return_value.reset(new WellNode(name, production_specification, injection_specification, phase_usage));
        } else {
            InjectionSpecification injection_specification;
            if (deck.hasField("GCONINJE")) {
                GCONINJE gconinje = deck.getGCONINJE();
                for (size_t i = 0; i < gconinje.gconinje.size(); i++) {
                    if (gconinje.gconinje[i].group_ == name) {
                        GconinjeLine line = gconinje.gconinje[i];
                        injection_specification.injector_type_ = toInjectorType(line.injector_type_);
                        injection_specification.control_mode_ = toInjectionControlMode(line.control_mode_);
                        injection_specification.surface_flow_max_rate_ = line.surface_flow_max_rate_;
                        injection_specification.reservoir_flow_max_rate_ = line.resv_flow_max_rate_;
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
                        production_specification.water_max_rate_ = line.water_max_rate_;
                        production_specification.gas_max_rate_ = line.gas_max_rate_;
                        production_specification.liquid_max_rate_ = line.liquid_max_rate_;
                        production_specification.procedure_ = toProductionProcedure(line.procedure_);
                    }
                }
            }

            return_value.reset(new WellsGroup(name, production_specification, injection_specification, phase_usage));
        }

        return return_value;
    }
}
