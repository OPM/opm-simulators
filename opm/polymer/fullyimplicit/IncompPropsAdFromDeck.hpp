/**/
#ifndef OPM_INCOMPPROPSADFROMDECK_HEADER_INCLUDED
#define OPM_INCOMPPROPSADFROMDECK_HEADER_INCLUDED
#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/core/props/rock/RockFromDeck.hpp>
#include <opm/core/props/pvt/PvtPropertiesIncompFromDeck.hpp>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

struct UnstructuredGrid;
namespace Opm
{
    class IncompPropsAdFromDeck : public IncompPropsAdInterface
    {
    public:
        IncompPropsAdFromDeck(const EclipseGridParser& deck,
                              const UnstructuredGrid&  grid);  
        ~IncompPropsAdFromDeck();

        //--Rock interface--
        int numDimensions() const;
        int numCells() const;
        const double* porosity() const;
        const double* permeability() const;

        // --Fluid interface--
        int numPhases() const;
        const double* viscosity() const;
        const double* density() const;
        const double* surfaceDensity() const;

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef std::vector<int> Cells;


        std::vector<V> relperm(const V& sw,
                               const V& so,
                               const Cells& cells) const;

        std::vector<ADB> relperm(const ADB& sw,
                                 const ADB& so,
                                 const Cells& cells) const;

    private:
        RockFromDeck rock_;
        PvtPropertiesIncompFromDeck pvt_;
        SaturationPropsFromDeck<SatFuncStone2Uniform> satprops_;
    };

} //namespace Opm

#endif// OPM_INCOMPPROPSADFROMDECK_HEADER_INCLUDED
