/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILPROPSADFROMDECK_HEADER_INCLUDED
#define OPM_BLACKOILPROPSADFROMDECK_HEADER_INCLUDED

#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/props/rock/RockFromDeck.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

namespace Opm
{

    class SinglePvtInterface;

    /// This class implements the AD-adapted fluid interface for
    /// three-phase black-oil. It requires an input deck from which it
    /// reads all relevant property data.
    ///
    /// Most methods are available in two overloaded versions, one
    /// taking a constant vector and returning the same, and one
    /// taking an AD type and returning the same. Derivatives are not
    /// returned separately by any method, only implicitly with the AD
    /// version of the methods.
    class BlackoilPropsAdFromDeck : public BlackoilPropsAdInterface
    {
    public:
        /// Constructor wrapping an opm-core black oil interface.
        BlackoilPropsAdFromDeck(const EclipseGridParser& deck,
                                const UnstructuredGrid& grid,
                                const bool init_rock = true );

        /// Constructor wrapping an opm-core black oil interface.
        BlackoilPropsAdFromDeck(Opm::DeckConstPtr newParserDeck,
                                const UnstructuredGrid& grid,
                                const bool init_rock = true );

        ////////////////////////////
        //      Rock interface    //
        ////////////////////////////

        /// \return   D, the number of spatial dimensions.
        int numDimensions() const;

        /// \return   N, the number of cells.
        int numCells() const;

        /// \return   Array of N porosity values.
        const double* porosity() const;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        const double* permeability() const;


        ////////////////////////////
        //      Fluid interface   //
        ////////////////////////////

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef std::vector<int> Cells;

        /// \return   Number of active phases (also the number of components).
        int numPhases() const;

        /// \return   Object describing the active phases.
        PhaseUsage phaseUsage() const;

        // ------ Canonical named indices for each phase ------

        /// Canonical named indices for each phase.
        enum PhaseIndex { Water = 0, Oil = 1, Gas = 2 };


        // ------ Density ------

        /// Densities of stock components at surface conditions.
        /// \return Array of 3 density values.
        const double* surfaceDensity() const;


        // ------ Viscosity ------

        /// Water viscosity.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        V muWat(const V& pw,
                const Cells& cells) const;

        /// Oil viscosity.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        V muOil(const V& po,
                const V& rs,
                const std::vector<PhasePresence>& cond,
                const Cells& cells) const;

        /// Gas viscosity.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        V muGas(const V& pg,
                const Cells& cells) const;

        /// Oil viscosity.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  rv     Array of n vapor oil/gas ratio
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        V muGas(const V& po,
                const V& rv,
                const std::vector<PhasePresence>& cond,
                const Cells& cells) const;

        /// Water viscosity.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muWat(const ADB& pw,
                  const Cells& cells) const;

        /// Oil viscosity.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muOil(const ADB& po,
                  const ADB& rs,
                  const std::vector<PhasePresence>& cond,
                  const Cells& cells) const;

        /// Gas viscosity.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muGas(const ADB& pg,
                  const Cells& cells) const;

        /// Gas viscosity.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muGas(const ADB& pg,
                  const ADB& rv,
                  const std::vector<PhasePresence>& cond,
                  const Cells& cells) const;

        // ------ Formation volume factor (b) ------

        /// Water formation volume factor.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        V bWat(const V& pw,
               const Cells& cells) const;

        /// Oil formation volume factor.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        V bOil(const V& po,
               const V& rs,
               const std::vector<PhasePresence>& cond,
               const Cells& cells) const;

        /// Gas formation volume factor.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        V bGas(const V& pg,
               const Cells& cells) const;

        /// Gas formation volume factor.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  rv     Array of n vapor oil/gas ratio
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        V bGas(const V& pg,
               const V& rv,
               const std::vector<PhasePresence>& cond,
               const Cells& cells) const;

        /// Water formation volume factor.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bWat(const ADB& pw,
                 const Cells& cells) const;

        /// Oil formation volume factor.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bOil(const ADB& po,
                 const ADB& rs,
                 const std::vector<PhasePresence>& cond,
                 const Cells& cells) const;

        /// Gas formation volume factor.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bGas(const ADB& pg,
                 const Cells& cells) const;

        /// Gas formation volume factor.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  rv     Array of n vapor oil/gas ratio
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bGas(const ADB& pg,
               const ADB& rv,
               const std::vector<PhasePresence>& cond,
               const Cells& cells) const;

        // ------ Rs bubble point curve ------

        /// Bubble point curve for Rs as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        V rsSat(const V& po,
                const Cells& cells) const;

        /// Bubble point curve for Rs as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        ADB rsSat(const ADB& po,
                  const Cells& cells) const;

        // ------ Rv condensation curve ------

        /// Condensation curve for Rv as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        V rvSat(const V& po,
                const Cells& cells) const;

        /// Condensation curve for Rv as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        ADB rvSat(const ADB& po,
                  const Cells& cells) const;

        // ------ Relative permeability ------

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  sg     Array of n gas saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 3 elements, each an array of n relperm values,
        ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
        std::vector<V> relperm(const V& sw,
                               const V& so,
                               const V& sg,
                               const Cells& cells) const;

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  sg     Array of n gas saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 3 elements, each an array of n relperm values,
        ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
        std::vector<ADB> relperm(const ADB& sw,
                                 const ADB& so,
                                 const ADB& sg,
                                 const Cells& cells) const;

        /// Capillary pressure for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  sg     Array of n gas saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 3 elements, each an array of n capillary pressure values,
        ///                    containing the offsets for each p_g, p_o, p_w. The capillary pressure between
        ///                    two arbitrary phases alpha and beta is then given as p_alpha - p_beta.
        std::vector<ADB> capPress(const ADB& sw,
                                  const ADB& so,
                                  const ADB& sg,
                                  const Cells& cells) const;
                                  
        /// Saturation update for hysteresis behavior.
        /// \param[in]  cells       Array of n cell indices to be associated with the saturation values.
        void updateSatHyst(const std::vector<double>& saturation,
                           const std::vector<int>& cells);

    private:
        RockFromDeck rock_;
        boost::scoped_ptr<SaturationPropsInterface> satprops_;
        PhaseUsage phase_usage_;
        std::vector<boost::shared_ptr<SinglePvtInterface> > props_;
        double densities_[BlackoilPhases::MaxNumPhases];
    };


} // namespace Opm

#endif // OPM_BLACKOILPROPSADFROMDECK_HEADER_INCLUDED
