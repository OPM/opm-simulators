#ifndef OPM_WELLSGROUP_HPP
#define	OPM_WELLSGROUP_HPP

#include <opm/core/InjectionSpecification.hpp>
#include <opm/core/ProductionSpecification.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <string>


namespace Opm
{
    // Need to forward declare this one, some of the methods in the base 
    // class returns pointers to it.
    class WellNode;
    
    /// Basic information needed for group control (each group should typically
    /// not exceed the sum of its leaf nodes)
    struct WellPhasesSummed {
        WellPhasesSummed();
        double bhp_sum;
        double rate_sum;
        
        /// Sums each component
        void operator+=(const WellPhasesSummed& other);
    };

    class WellsGroupInterface
    {
    public:
        WellsGroupInterface(const std::string& name,
                            ProductionSpecification prod_spec,
                            InjectionSpecification inj_spec);
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
        
        
        /// \returns true if the object is a leaf node (WellNode), false otherwise.
        virtual bool isLeafNode() const;
        
        /// \returns the pointer to the WellsGroupInterface with the given name. NULL if 
        ///          the name is not found.a
        virtual WellsGroupInterface* findGroup(std::string name_of_node) = 0;

        /// Sets the parent
        /// \param[in] parent the pointer to the parent
        void setParent(WellsGroupInterface* parent);
        
        /// Gets the parent of the group, NULL if no parent.
        const WellsGroupInterface* getParent() const;
        
        /// Recursively calculate the guide rate for each member of the well group.
        /// This should be called after the guide rates are set to the non-normalized values.
        virtual void calculateGuideRates() = 0;
        
        /// Calculates the number of leaf nodes in the given group. 
        /// A leaf node is defined to have one leaf node in its group.
        virtual int numberOfLeafNodes() = 0;

        /// Checks if each condition is met, applies well controls where needed
        /// (that is, it either changes the active control of violating wells, or shuts
        /// down wells). Only one change is applied per invocation. Typical use will be
        /// \code
        /// solve_pressure();
        /// while(!group.conditionsMet(well_bhp, well_rate, summed_phases)) {
        ///     solve_pressure();
        /// }
        /// \endcode
        ///
        /// \note It's highly recommended to use the conditionsMet found in WellsManager.
        /// \param[in]    well_bhp  A vector containing the bhp for each well. Is assumed 
        ///                         to be ordered the same way as the related Wells-struct.
        /// \param[in]    well_rate A vector containing the rate for each well. Is assumed 
        ///                         to be ordered the same way as the related Wells-struct.
        /// \param[out]   summed_phases Will at end of invocation contain the summed phases
        ///                             (bhp, rate ,etc.) for the group.
        /// \param[in]    epsilon   The error tolerated for each inequality. Formally, it will accept
        ///                         (a - b <= epsilon) as (a <= b).
        /// \return true if no violations were found, false otherwise (false also implies a change).
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_rate,
                                   WellPhasesSummed& summed_phases,
                                   const double epsilon = 1e-8) = 0;
        
        /// Sets the current active control to the provided one for all wells within the group
        /// \note Also changes the target based on type.
        /// \param[in] type the type to change to which the control is changed.
        virtual void applyControl(const WellControlType type) = 0;
        
        /// Gets the worst offending well based on the input
        /// \param values A vector of a values for each well. This is assumed to be ordered the same way as the 
        ///               relevant Wells struct.
        /// \return first will be a pointer to the worst offending well, second will be the obtained value at that well.
        virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& values) = 0;
        
    protected:
           WellsGroupInterface* parent_;

    private:
        std::string name_;
        ProductionSpecification production_specification_;
        InjectionSpecification injection_specification_;
    };



    class WellsGroup : public WellsGroupInterface
    {
    public:
        WellsGroup(const std::string& name,
                   ProductionSpecification prod_spec,
                   InjectionSpecification inj_spec);

        virtual WellsGroupInterface* findGroup(std::string name_of_node);

        void addChild(std::tr1::shared_ptr<WellsGroupInterface> child);
        
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_rate,
                                   WellPhasesSummed& summed_phases,
                                   const double epsilon = 1e-8);
        
        
        virtual void calculateGuideRates();

        virtual int numberOfLeafNodes();
        virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& values);
        virtual void applyControl(const WellControlType type);

    private:
        std::vector<std::tr1::shared_ptr<WellsGroupInterface> > children_;
    };



    class WellNode : public WellsGroupInterface
    {
    public:
        WellNode(const std::string& name,
                 ProductionSpecification prod_spec,
                InjectionSpecification inj_spec);

        virtual WellsGroupInterface* findGroup(std::string name_of_node);
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_rate,
                                   WellPhasesSummed& summed_phases,
                                   const double epsilon = 1e-8);
        
        virtual bool isLeafNode() const;
        
        void setWellsPointer(Wells* wells, int self_index);
        
        virtual void calculateGuideRates();
        virtual int numberOfLeafNodes();
        
        // Shuts the well (in the well struct)
        void shutWell();
        
        virtual std::pair<WellNode*, double> getWorstOffending(const std::vector<double>& values);
        virtual void applyControl(const WellControlType type);
        
    private:
        Wells* wells_;
        int self_index_;
    };

    /// Doc me!
    std::tr1::shared_ptr<WellsGroupInterface> createWellsGroup(std::string name, 
                                                               const EclipseGridParser& deck);


}
#endif	/* OPM_WELLSGROUP_HPP */

