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

#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/fluid/blackoil/BlackoilPhases.hpp>
#include <vector>

namespace Opm
{

    class SaturationPropsFromDeck : public BlackoilPhases
    {
    public:
        /// Default constructor.
        SaturationPropsFromDeck();

        /// Initialize from deck.
        /// global_cell maps from grid cells to their original logical Cartesian indices.
        void init(const EclipseGridParser& deck,
                  const std::vector<int>& global_cell);

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
        class SatFuncSet
        {
        public:
            void init(const EclipseGridParser& deck, const int table_num, PhaseUsage phase_usg);
            void evalKr(const double* s, double* kr) const;
            void evalKrDeriv(const double* s, double* kr, double* dkrds) const;
            void evalPc(const double* s, double* pc) const;
            void evalPcDeriv(const double* s, double* pc, double* dpcds) const;
            double smin_[PhaseUsage::MaxNumPhases];
            double smax_[PhaseUsage::MaxNumPhases];
        private:
            PhaseUsage phase_usage; // A copy of the outer class' phase_usage_.
            UniformTableLinear<double> krw_;
            UniformTableLinear<double> krow_;
            UniformTableLinear<double> pcow_;
            UniformTableLinear<double> krg_;
            UniformTableLinear<double> krog_;
            UniformTableLinear<double> pcog_;
            double krocw_; // = krow_(s_wc)
        };
        std::vector<SatFuncSet> satfuncset_;
        std::vector<int> cell_to_func_; // = SATNUM - 1

        const SatFuncSet& funcForCell(const int cell) const;
    };



} // namespace Opm




#endif // OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
