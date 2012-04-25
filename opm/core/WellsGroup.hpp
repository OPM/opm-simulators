#ifndef OPM_WELLSGROUP_HPP
#define	OPM_WELLSGROUP_HPP

#include <opm/core/InjectionSpecification.hpp>
#include <opm/core/ProductionSpecification.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/grid.h>
#include <string>


namespace Opm
{
    class WellNode;
    struct WellPhasesSummed {
        WellPhasesSummed();
        double bhp_sum;
        double rate_sum;
        
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
        
        void setParent(WellsGroupInterface* parent);
        const WellsGroupInterface* getParent() const;
        
        /// Recursively calculate the guide rate for each member of the well group.
        /// This should be called after the guide rates are set to the non-normalized values.
        virtual void calculateGuideRates() = 0;
        
        /// Calculates the number of leaf nodes in the given group. 
        /// A leaf node is defined to have one leaf node in its group.
        virtual int numberOfLeafNodes() = 0;

        /// Fills the WellControlResult parameter with all exceed information
        virtual bool conditionsMet(const std::vector<double>& well_bhp,
                                   const std::vector<double>& well_rate,
                                   WellPhasesSummed& summed_phases,
                                   const double epsilon = 1e-8) = 0;
        
        virtual void applyControl(const WellControlType type) = 0;
        
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

