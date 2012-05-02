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

#ifndef OPM_WELLSGROUP_HPP
#define	OPM_WELLSGROUP_HPP

#include <opm/core/InjectionSpecification.hpp>
#include <opm/core/ProductionSpecification.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <opm/core/fluid/blackoil/BlackoilPhases.hpp>
#include <string>


namespace Opm
{
    // Need to forward declare this one, some of the methods in the base 
    // class returns pointers to it.
    class WellNode;
    
    /// Basic information needed for group control (each group should typically
    /// not exceed the sum of its leaf nodes)
    struct WellPhasesSummed
    {
        WellPhasesSummed();
        double res_inj_rates[3];
        double res_prod_rates[3];
        double surf_inj_rates[3];
        double surf_prod_rates[3];

        /// Sums each component
        void operator+=(const WellPhasesSummed& other);
    };

    class WellsGroupInterface
    {
    public:
        WellsGroupInterface(const std::string& name,
                            const ProductionSpecification& prod_spec,
                            const InjectionSpecification& inj_spec,
                            const PhaseUsage& phase_usage);
        virtual ~WellsGroupInterface();

        /// The unique identifier for the well or well group.
        const std::string& name();
        
        /// Production specifications for the well or well group.
        const ProductionSpecification& prodSpec() const;
        
        /// Injection specifications for the well or well group.
        const InjectionSpecification& injSpec() const;
        
        /// Production specifications for the well or well group.
        ProductionSpecification& prodSpec();
        
        /// Injection specifications for the well or well group.
        InjectionSpecification& injSpec();

        /// Phase usage information.
        const PhaseUsage& phaseUsage() const;
        
        /// \returns true if the object is a leaf node (WellNode), false otherwise.
        virtual bool isLeafNode() const;
        
        /// \returns the pointer to the WellsGroupInterface with the given name. NULL if 
        ///          the name is not found.a
        virtual WellsGroupInterface* findGroup(const std::string& name_of_node) = 0;

        /// Sets the parent
        /// \param[in] parent the pointer to the parent
        void setParent(WellsGroupInterface* parent);
        
        /// Gets the parent of the group, NULL if no parent.
        const WellsGroupInterface* getParent() const;
        
        /// Calculates the number of leaf nodes in the given group. 
        /// A leaf node is defined to have one leaf node in its group.
        virtual int numberOfLeafNodes() = 0;

        /// Checks if each condition is met, applies well controls where needed
        /// (that is, it either changes the active control of violating wells, or shuts
        /// down wells). Only one change is applied per invocation. Typical use will be
        /// \code
        /// solve_pressure();
        /// while(!group.conditionsMet(...)) {
        ///     solve_pressure();
        /// }
        /// \endcode
        ///
        /// \note It's highly recommended to use the conditionsMet found in WellsManager.
        /// \param[in]    well_bhp  A vector containing the bhp for each well. Is assumed 
        ///                         to be ordered the same way as the related Wells-struct.
        /// \param[in]    well_reservoirrates_phase
        ///                         A vector containing reservoir rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[in]    well_surfacerates_phase
        ///                         A vector containing surface rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[out]   summed_phases Will at end of invocation contain the summed phase rates
        ///                             (rate ,etc.) for the group.
        /// \return true if no violations were found, false otherwise (false also implies a change).
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_reservoirrates_phase,
                                   const std::vector<double>& well_surfacerates_phase,
                                   WellPhasesSummed& summed_phases) = 0;
        
        /// Sets the current active control to the provided one for all injectors within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
        /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyInjGroupControl(const InjectionSpecification::ControlMode control_mode,
                                          const double target,
                                          const bool forced) = 0;
        /// Sets the current active control to the provided one for all producers within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
        /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyProdGroupControl(const ProductionSpecification::ControlMode control_mode,
                                           const double target,
                                           const bool forced) = 0;

        /// Gets the worst offending well based on the input
        /// \param[in]    well_reservoirrates_phase
        ///                         A vector containing reservoir rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[in]    well_surfacerates_phase
        ///                         A vector containing surface rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[in]   mode
        ///                         The relevant control mode to find the maximum over.
        /// \return first will be a pointer to the worst offending well, second will be the obtained value at that well.
        virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& well_reservoirrates_phase,
                                                               const std::vector<double>& well_surfacerates_phase,
                                                               ProductionSpecification::ControlMode mode) = 0;
        
        /// Gets the target rate for the given mode.
        double getTarget(ProductionSpecification::ControlMode mode);
        
        /// Gets the target rate for the given mode.
        double getTarget(InjectionSpecification::ControlMode mode);
        
        /// Applies any production group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyProdGroupControls() = 0;
        
        /// Applies any injection group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyInjGroupControls() = 0;
        
        /// Calculates the production guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double productionGuideRate(bool only_group) = 0;
        
        /// Calculates the injection guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double injectionGuideRate(bool only_group) = 0;
        
    protected:
        /// Calculates the correct rate for the given ProductionSpecification::ControlMode
        double rateByMode(const double* res_rates, 
                          const double* surf_rates,
                          const ProductionSpecification::ControlMode mode);

        /// Calculates the correct rate for the given InjectionSpecification::ControlMode
        double rateByMode(const double* res_rates, 
                          const double* surf_rates,
                          const InjectionSpecification::ControlMode mode);

        WellsGroupInterface* parent_;

    private:
        std::string name_;
        ProductionSpecification production_specification_;
        InjectionSpecification injection_specification_;
        PhaseUsage phase_usage_;
    };



    class WellsGroup : public WellsGroupInterface
    {
    public:
        WellsGroup(const std::string& name,
                   const ProductionSpecification& prod_spec,
                   const InjectionSpecification& inj_spec,
                   const PhaseUsage& phase_usage);

        virtual WellsGroupInterface* findGroup(const std::string& name_of_node);

        void addChild(std::tr1::shared_ptr<WellsGroupInterface> child);
        
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_reservoirrates_phase,
                                   const std::vector<double>& well_surfacerates_phase,
                                   WellPhasesSummed& summed_phases);
        
        virtual int numberOfLeafNodes();
        virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& well_reservoirrates_phase,
                                                               const std::vector<double>& well_surfacerates_phase,
                                                               ProductionSpecification::ControlMode mode);

        /// Sets the current active control to the provided one for all injectors within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
        /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyInjGroupControl(const InjectionSpecification::ControlMode control_mode,
                                          const double target,
                                          bool forced);

        /// Sets the current active control to the provided one for all producers within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
        /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyProdGroupControl(const ProductionSpecification::ControlMode control_mode,
                                           const double target,
                                           bool forced);
        
        /// Applies any production group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyProdGroupControls();
        
        /// Applies any injection group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyInjGroupControls();
        
        /// Calculates the production guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double productionGuideRate(bool only_group);
        
        /// Calculates the injection guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double injectionGuideRate(bool only_group);

    private:
        std::vector<std::tr1::shared_ptr<WellsGroupInterface> > children_;
    };



    class WellNode : public WellsGroupInterface
    {
    public:
        WellNode(const std::string& name,
                 const ProductionSpecification& prod_spec,
                 const InjectionSpecification& inj_spec,
                 const PhaseUsage& phase_usage);

        virtual WellsGroupInterface* findGroup(const std::string& name_of_node);
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_reservoirrates_phase,
                                   const std::vector<double>& well_surfacerates_phase,
                                   WellPhasesSummed& summed_phases);
        
        virtual bool isLeafNode() const;
        
        void setWellsPointer(Wells* wells, int self_index);
        
        virtual int numberOfLeafNodes();
        
        // Shuts the well (in the well struct)
        void shutWell();
        
       virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& well_reservoirrates_phase,
                                                               const std::vector<double>& well_surfacerates_phase,
                                                               ProductionSpecification::ControlMode mode);

        /// Sets the current active control to the provided one for all injectors within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
       /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyInjGroupControl(const InjectionSpecification::ControlMode control_mode,
                                          const double target,
                                          bool forced);

        /// Sets the current active control to the provided one for all producers within the group.
        /// After this call, the combined rate (which rate depending on control_mode) of the group
        /// shall be equal to target.
        /// \param[in] forced if true, all children will be set under group control, otherwise
        ///                   only children that are under group control will be changed.
        virtual void applyProdGroupControl(const ProductionSpecification::ControlMode control_mode,
                                           const double target,
                                           bool forced);
        
        /// Applies any production group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyProdGroupControls();
        
        /// Applies any injection group control relevant to all children nodes.
        /// If no group control is set, this is called recursively to the children.
        virtual void applyInjGroupControls();
        
        /// Calculates the production guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double productionGuideRate(bool only_group);
        
        /// Calculates the injection guide rate for the group.
        /// \param[in] only_group If true, will only accumelate guide rates for 
        ///                       wells under group control
        virtual double injectionGuideRate(bool only_group);

    private:
        Wells* wells_;
        int self_index_;
        int group_control_index_;
        bool shut_well_;
    };

    /// Creates the WellsGroupInterface for the given name
    /// \param[in] name the name of the wells group.
    /// \param[in] deck the deck from which to fetch information.
    std::tr1::shared_ptr<WellsGroupInterface> createWellsGroup(const std::string& name, 
                                                               const EclipseGridParser& deck);


}
#endif	/* OPM_WELLSGROUP_HPP */

