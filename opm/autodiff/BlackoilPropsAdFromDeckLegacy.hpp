/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU.

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

#ifndef OPM_BLACKOILPROPSADFROMDECKLEGACY_HEADER_INCLUDED
#define OPM_BLACKOILPROPSADFROMDECKLEGACY_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>

namespace Opm
{
    /// This class implements the AD-adapted fluid interface for
    /// three-phase black-oil. It requires an input deck from which it
    /// reads all relevant property data.
    ///
    /// Most methods are available in two overloaded versions, one
    /// taking a constant vector and returning the same, and one
    /// taking an AD type and returning the same. Derivatives are not
    /// returned separately by any method, only implicitly with the AD
    /// version of the methods.
    class BlackoilPropsAdFromDeckLegacy : public BlackoilPropsAdFromDeck
    {
    public:
        BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                      const Opm::EclipseState& eclState,
                                      std::shared_ptr<MaterialLawManager> materialLawManager,
                                      const UnstructuredGrid& grid,
                                      const bool init_rock = true );

#ifdef HAVE_OPM_GRID
        /// Constructor to create a blackoil properties from an ECL deck.
        ///
        /// The materialLawManager parameter represents the object from opm-material
        /// which handles the creating and updating parameter objects for the capillary
        /// pressure/relperm relations for each grid cell. This object is created
        /// internally for the constructors below, but if it is already available
        /// externally some performance can be gained by creating it only once.
        ///
        /// \param deck The unprocessed ECL deck from opm-parser
        /// \param eclState The processed ECL deck from opm-parser
        /// \param materialLawManager The container for the material law parameter objects
        /// \param grid The grid upon which the simulation is run on.
        /// \param init_rock If true the rock properties (rock compressibility and
        ///                  reference pressure) are read from the deck
        BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                      const Opm::EclipseState& eclState,
                                      std::shared_ptr<MaterialLawManager> materialLawManager,
                                      const Dune::CpGrid& grid,
                                      const bool init_rock = true );
#endif

        /// Constructor to create a blackoil properties from an ECL deck.
        ///
        /// \param deck The unprocessed ECL deck from opm-parser
        /// \param eclState The processed ECL deck from opm-parser
        /// \param grid The grid upon which the simulation is run on.
        /// \param init_rock If true the rock properties (rock compressibility and
        ///                  reference pressure) are read from the deck
        BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                      const Opm::EclipseState& eclState,
                                      const UnstructuredGrid& grid,
                                      const bool init_rock = true );

#ifdef HAVE_OPM_GRID
        /// Constructor to create a blackoil properties from an ECL deck.
        ///
        /// \param deck The unprocessed ECL deck from opm-parser
        /// \param eclState The processed ECL deck from opm-parser
        /// \param grid The grid upon which the simulation is run on.
        /// \param init_rock If true the rock properties (rock compressibility and
        ///                  reference pressure) are read from the deck
        BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                      const Opm::EclipseState& eclState,
                                      const Dune::CpGrid& grid,
                                      const bool init_rock = true );
#endif

        /// \brief Constructor to create properties for a subgrid
        ///
        /// This copies all properties that are not dependant on the
        /// grid size from an existing properties object
        /// and the number of cells. All properties that do not depend
        /// on the grid dimension will be copied. For the rest will have
        /// the correct size but the values will be undefined.
        ///
        /// \param props            The property object to copy from.
        /// \param materialLawManager The container for the material law parameter objects.
        ///                           Initialized only for the subgrid
        /// \param number_of_cells  The number of cells of the subgrid.
        BlackoilPropsAdFromDeckLegacy(const BlackoilPropsAdFromDeckLegacy& props,
                                std::shared_ptr<MaterialLawManager> materialLawManager,
                                const int number_of_cells);

        ////////////////////////////
        //      Fluid interface   //
        ////////////////////////////

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;

        // ------ Density ------

        /// \return Array of n density values for phase given by phaseIdx.
        V surfaceDensity(const int phaseIdx, const Cells& cells) const;

        // ------ Viscosity ------

        /// Water viscosity.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muWat(const ADB& pw,
                  const ADB& T,
                  const Cells& cells) const;

        /// Oil viscosity.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muOil(const ADB& po,
                  const ADB& T,
                  const ADB& rs,
                  const std::vector<PhasePresence>& cond,
                  const Cells& cells) const;

        /// Gas viscosity.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  rv     Array of n vapor oil/gas ratios.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n viscosity values.
        ADB muGas(const ADB& pg,
                  const ADB& T,
                  const ADB& rv,
                  const std::vector<PhasePresence>& cond,
                  const Cells& cells) const;

        // ------ Formation volume factor (b) ------

        /// Water formation volume factor.
        /// \param[in]  pw     Array of n water pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bWat(const ADB& pw,
                 const ADB& T,
                 const Cells& cells) const;

        /// Oil formation volume factor.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  rs     Array of n gas solution factor values.
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bOil(const ADB& po,
                 const ADB& T,
                 const ADB& rs,
                 const std::vector<PhasePresence>& cond,
                 const Cells& cells) const;

        /// Gas formation volume factor.
        /// \param[in]  pg     Array of n gas pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  rv     Array of n vapor oil/gas ratio
        /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n formation volume factor values.
        ADB bGas(const ADB& pg,
                 const ADB& T,
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
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        V rsSat(const V& po,
                const V& so,
                const Cells& cells) const;

        /// Bubble point curve for Rs as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        ADB rsSat(const ADB& po,
                  const Cells& cells) const;

        /// Bubble point curve for Rs as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n bubble point values for Rs.
        ADB rsSat(const ADB& po,
                  const ADB& so,
                  const Cells& cells) const;

        // ------ Rv condensation curve ------

        /// Condensation curve for Rv as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n condensation point values for Rv.
        ADB rvSat(const ADB& po,
                  const Cells& cells) const;

        /// Condensation curve for Rv as function of oil pressure.
        /// \param[in]  po     Array of n oil pressure values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
        /// \return            Array of n condensation point values for Rv.
        ADB rvSat(const ADB& po,
                  const ADB& so,
                  const Cells& cells) const;

        // ------ Relative permeability ------

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
                                  

        /// Returns the bubble point pressures
        std::vector<double> bubblePointPressure(const Cells& cells,
                const V& T,
                const V& rs) const;

        /// Returns the dew point pressures
        std::vector<double> dewPointPressure(const Cells& cells,
                const V& T,
                const V& rv) const;

        /// Obtain the scaled critical oil in gas saturation values.
        /// \param[in]  cells  Array of cell indices.
        /// \return Array of critical oil in gas saturaion values.
        V scaledCriticalOilinGasSaturations(const Cells& cells) const;

        V scaledCriticalGasSaturations(const Cells& cells) const;


    private:
        /// Correction to rs/rv according to kw VAPPARS
        void applyVap(V& r,
                      const V& so,
                      const std::vector<int>& cells,
                      const double vap) const;

        void applyVap(ADB& r,
                      const ADB& so,
                      const std::vector<int>& cells,
                      const double vap) const;
    };
} // namespace Opm

#endif // OPM_BLACKOILPROPSADFROMDECKLEGACY_HEADER_INCLUDED
