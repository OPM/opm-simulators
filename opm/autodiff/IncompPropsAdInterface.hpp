/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#ifndef OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED
#define OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>

namespace Opm
{
    class IncompPropsAdInterface
    {
    public:
        /// Virtual destructor for inheritance.
        virtual ~IncompPropsAdInterface();

        ////////////////////////////
        //      Rock interface    //
        ////////////////////////////

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const = 0;

        /// \return   N, the number of cells.
        virtual int numCells() const = 0;

        /// \return   Array of N porosity values.
        virtual const double* porosity() const = 0;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const = 0;

        ////////////////////////////
        //      Fluid interface   //
        ////////////////////////////

        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef std::vector<int> Cells;

        /// \return   P, Number of active phases (also the number of components).
        virtual int numPhases() const = 0;

        // ------ Density ------

        /// Densities of stock components at surface conditions.
        /// \return Array of P density values.
        virtual const double* surfaceDensity() const = 0;

        // ------ Viscosity ------

        /// Viscosity of stock components at surface conditions.
        /// \return Array of P viscosity values.
        virtual const double* viscosity() const = 0;

        // ------ Relative permeability ------

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro.
        virtual
        std::vector<V> relperm(const V& sw,
                               const V& so,
                               const Cells& cells) const = 0;

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro.
        virtual
        std::vector<ADB> relperm(const ADB& sw,
                                 const ADB& so,
                                 const Cells& cells) const = 0;

        /// Capillary pressure for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n capillary pressure values,
        ///                    containing the offsets for each p_o, p_w. The capillary pressure between
        ///                    two arbitrary phases alpha and beta is then given as p_alpha - p_beta.
        virtual
        std::vector<ADB> capPress(const ADB& sw,
                                  const ADB& so,
                                  const Cells& cells) const = 0;

    };
} // namespace Opm

#endif// OPM_INCOMPPROPSADINTERFACE_HEADER_INCLUDED
