/**/

#ifndef OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED
#define OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED

#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/core/props/rock/RockBasic.hpp>
#include <opm/core/props/pvt/PvtPropertiesBasic.hpp>
#include <opm/core/props/satfunc/SaturationPropsBasic.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
namespace Opm
{
    class IncompPropsAdBasic : public IncompPropsAdInterface
    {
    public:
        IncompPropsAdBasic(const parameter::ParameterGroup& param,
                           const int dim,
                           const int num_cells);
        IncompPropsAdBasic(const int num_phases,
                              const SaturationPropsBasic::RelPermFunc& relpermfunc,
                              const std::vector<double>& rho,
                              const std::vector<double>& mu,
                              const double porosity,
                              const double permeability,
                              const int dim,
                              const int num_cells);

         ~IncompPropsAdBasic();
         int numDimensions() const;
         int numCells() const;
         const double* porosity() const;
         const double* permeability() const;
         
         int numPhases() const;
         const double* viscosity() const;
         const double* density() const;
         const double* surfaceDensity() const;

         typedef AutoDiffBlock<double> ADB;
         typedef ADB::V V;
         std::vector<V> relperm(const V& sw,
                                const V& so,
                                const std::vector<int>& cells) const;
         std::vector<ADB> relperm(const ADB& sw,
                     const ADB& so,
                     const std::vector<int>& cells) const;
    private:
        RockBasic rock_;
        PvtPropertiesBasic pvt_;
        SaturationPropsBasic satprops_;
        std::vector<double> viscosity_;
    };
}
#endif // OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED
