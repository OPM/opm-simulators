/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL

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

#ifndef OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED
#define OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED

#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/core/props/rock/RockBasic.hpp>
#include <opm/core/props/pvt/PvtPropertiesBasic.hpp>
#include <opm/core/props/satfunc/SaturationPropsBasic.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>

namespace Opm
{
    /// This class implements the AD-adapted fluid interface for
    /// two-phase oil-water. It requires an input deck from which it
    /// reads all relevant property data.
    ///
    /// Most methods are available in two overloaded versions, one
    /// taking a constant vector and returning the same, and one
    /// taking an AD type and returning the same. Derivatives are not
    /// returned separately by any method, only implicitly with the AD
    /// version of the methods.
    class IncompPropsAdBasic : public IncompPropsAdInterface
    {
    public:
		/// Constructor.
        IncompPropsAdBasic(const parameter::ParameterGroup& param,
                           const int dim,
                           const int num_cells);

		/// Constructor.
        IncompPropsAdBasic(const int num_phases,
                           const SaturationPropsBasic::RelPermFunc& relpermfunc,
                           const std::vector<double>& rho,
                           const std::vector<double>& mu,
                           const double porosity,
                           const double permeability,
                           const int dim,
                           const int num_cells);
 		/// Destructor.
        ~IncompPropsAdBasic();

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

        /// \return   P, Number of active phases (also the number of components).
        int numPhases() const;

        /// \return   Array of P viscosity values.
        const double* viscosity() const;

		/// Densities of fluid phases at reservoir conditions.
        /// \return   Array of P density values.
        const double* density() const;

		/// Densities of fluid phases at surface conditions.
        /// \return   Array of P density values.
        const double* surfaceDensity() const;

        // ------ Relative permeability ------
        
        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro.
        std::vector<V> relperm(const V& sw,
                               const V& so,
                               const std::vector<int>& cells) const;

        /// Relative permeabilities for all phases.
        /// \param[in]  sw     Array of n water saturation values.
        /// \param[in]  so     Array of n oil saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
        /// \return            An std::vector with 2 elements, each an array of n relperm values,
        ///                    containing krw, kro.
        std::vector<ADB> relperm(const ADB& sw,
        						 const ADB& so,
                     			 const std::vector<int>& cells) const;
    private:
        RockBasic rock_;
        PvtPropertiesBasic pvt_;
        SaturationPropsBasic satprops_;
        std::vector<double> viscosity_;
    };
} //namespace Opm

#endif // OPM_INCOMPPROPSADBASIC_HEADER_INCLUDED
