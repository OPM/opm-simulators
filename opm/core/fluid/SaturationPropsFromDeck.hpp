/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
#define OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/fluid/blackoil/BlackoilPhases.hpp>
#include <opm/core/fluid/SatFuncStone2.hpp>
#include <opm/core/fluid/SatFuncSimple.hpp>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class SaturationPropsFromDeck : public BlackoilPhases
    {
    public:
        /// Default constructor.
        SaturationPropsFromDeck();

        /// Initialize from deck and grid.
        /// \param  deck         Deck input parser
        /// \param  grid         Grid to which property object applies, needed for the 
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        void init(const EclipseGridParser& deck,
                  const UnstructuredGrid& grid);

        void init(const EclipseGridParser& deck,
                  const UnstructuredGrid& grid,
                  const parameter::ParameterGroup& param);

        /// \return   P, the number of phases.
        int numPhases() const;

        /// Relative permeability.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        void relperm(const int n,
                     const double* s,
                     const int* cells,
                     double* kr,
                     double* dkrds) const;

        /// Capillary pressure.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dpc_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        void capPress(const int n,
                      const double* s,
                      const int* cells,
                      double* pc,
                      double* dpcds) const;

	/// Obtain the range of allowable saturation values.
        /// \param[in]  n      Number of data points.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
	void satRange(const int n,
                      const int* cells,
		      double* smin,
		      double* smax) const;

    private:
        PhaseUsage phase_usage_;
        typedef SatFuncSimple satfunc_t;
        std::vector<satfunc_t> satfuncset_;
        std::vector<int> cell_to_func_; // = SATNUM - 1

        const satfunc_t& funcForCell(const int cell) const;
    };



} // namespace Opm




#endif // OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
