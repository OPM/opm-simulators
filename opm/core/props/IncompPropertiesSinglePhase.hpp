/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_INCOMPPROPERTIESSINGLEPHASE_HEADER_INCLUDED
#define OPM_INCOMPPROPERTIESSINGLEPHASE_HEADER_INCLUDED



#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/props/rock/RockFromDeck.hpp>

struct UnstructuredGrid;

namespace Opm
{

    /// Concrete class implementing the incompressible property
    /// interface for a simplified single-phase setting, reading all
    /// data and properties from eclipse deck input. The oil phase
    /// properties are used where applicable and available.
    ///
    /// Supports variable number of spatial dimensions, called D.
    /// Supports a single phase only.
    /// In general, when arguments call for n values of some vector or
    /// matrix property, such as saturation, they shall always be
    /// ordered cellwise:
    ///   [s^1_0 s^2_0 s^3_0 s^1_1 s^2_2 ... ]
    /// in which s^i_j denotes saturation of phase i in cell j.
    class IncompPropertiesSinglePhase : public IncompPropertiesInterface
    {
    public:
        /// Initialize from deck and grid.
        /// \param  deck         Deck input parser
        /// \param  eclState        The EclipseState (processed deck) produced by the opm-parser code
        /// \param  grid         Grid to which property object applies, needed for the
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        IncompPropertiesSinglePhase(Opm::DeckConstPtr deck,
                                 Opm::EclipseStateConstPtr eclState,
                                 const UnstructuredGrid& grid);

        /// Destructor.
        virtual ~IncompPropertiesSinglePhase();

        // ---- Rock interface ----

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const;

        /// \return   N, the number of cells.
        virtual int numCells() const;

        /// \return   Array of N porosity values.
        virtual const double* porosity() const;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const;


        // ---- Fluid interface ----

        /// \return   P, the number of phases (= 1).
        virtual int numPhases() const;

        /// \return Array of P (= 1) viscosity values.
        virtual const double* viscosity() const;

        /// Densities of fluid at reservoir conditions.
        /// \return Array of P (= 1) density values.
        virtual const double* density() const;

        /// Densities of fluid phases at surface conditions.
        /// \return Array of P (= 1) density values.
        virtual const double* surfaceDensity() const;

        /// Relative permeability. Always returns 1 (and 0 for derivatives).
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of n saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] kr     Array of n relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of n relperm derivative values,
        ///                    array must be valid before calling.
        virtual void relperm(const int n,
                             const double* s,
                             const int* cells,
                             double* kr,
                             double* dkrds) const;

        /// Capillary pressure. Always returns zero.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of n saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] pc     Array of n capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of n derivative values,
        ///                    array must be valid before calling.
        virtual void capPress(const int n,
                              const double* s,
                              const int* cells,
                              double* pc,
                              double* dpcds) const;


        /// Obtain the range of allowable saturation values.
        /// Saturation range is just the point 1 for this class
        /// \param[in]  n      Number of data points.
        /// \param[in]  cells  Array of n cell indices.
        /// \param[out] smin   Array of n minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of n maximum s values, array must be valid before calling.
        virtual void satRange(const int n,
                              const int* cells,
                              double* smin,
                              double* smax) const;
    private:
        RockFromDeck rock_;
        double surface_density_;
        double reservoir_density_;
        double viscosity_;
    };



} // namespace Opm



#endif // OPM_INCOMPPROPERTIESSINGLEPHASE_HEADER_INCLUDED
