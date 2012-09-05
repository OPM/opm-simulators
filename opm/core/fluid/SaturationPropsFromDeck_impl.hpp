/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
#define OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED


#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/utility/NonuniformTableLinear.hpp>
#include <opm/core/fluid/blackoil/phaseUsageFromDeck.hpp>
#include <opm/core/grid.h>

namespace Opm
{


    // ----------- Methods of SaturationPropsFromDeck ---------


    /// Default constructor.
    template <class SatFuncSet>
    SaturationPropsFromDeck<SatFuncSet>::SaturationPropsFromDeck()
    {
    }

    /// Initialize from deck.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::init(const EclipseGridParser& deck,
                                                   const UnstructuredGrid& grid,
                                                   const int samples)
    {
        phase_usage_ = phaseUsageFromDeck(deck);

        // Extract input data.
        // Oil phase should be active.
        if (!phase_usage_.phase_used[Liquid]) {
            THROW("SaturationPropsFromDeck::init()   --  oil phase must be active.");
        }

        // Obtain SATNUM, if it exists, and create cell_to_func_.
        // Otherwise, let the cell_to_func_ mapping be just empty.
        int satfuncs_expected = 1;
        if (deck.hasField("SATNUM")) {
            const std::vector<int>& satnum = deck.getIntegerValue("SATNUM");
            satfuncs_expected = *std::max_element(satnum.begin(), satnum.end());
            const int num_cells = grid.number_of_cells;
            cell_to_func_.resize(num_cells);
            const int* gc = grid.global_cell;
            for (int cell = 0; cell < num_cells; ++cell) {
                const int deck_pos = (gc == NULL) ? cell : gc[cell];
                cell_to_func_[cell] = satnum[deck_pos] - 1;
            }
        }

        // Find number of tables, check for consistency.
        enum { Uninitialized = -1 };
        int num_tables = Uninitialized;
        if (phase_usage_.phase_used[Aqua]) {
            const SWOF::table_t& swof_table = deck.getSWOF().swof_;
            num_tables = swof_table.size();
            if (num_tables < satfuncs_expected) {
                THROW("Found " << num_tables << " SWOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
        }
        if (phase_usage_.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            int num_sgof_tables = sgof_table.size();
            if (num_sgof_tables < satfuncs_expected) {
                THROW("Found " << num_tables << " SGOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
            if (num_tables == Uninitialized) {
                num_tables = num_sgof_tables;
            } else if (num_tables != num_sgof_tables) {
                THROW("Inconsistent number of tables in SWOF and SGOF.");
            }
        }

        // Initialize tables.
        satfuncset_.resize(num_tables);
        for (int table = 0; table < num_tables; ++table) {
            satfuncset_[table].init(deck, table, phase_usage_, samples);
        }
    }




    /// \return   P, the number of phases.
    template <class SatFuncSet>
    int SaturationPropsFromDeck<SatFuncSet>::numPhases() const
    {
        return phase_usage_.num_phases;
    }




    /// Relative permeability.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::relperm(const int n,
                                          const double* s,
                                          const int* cells,
                                          double* kr,
                                          double* dkrds) const
    {
        ASSERT (cells != 0);

        const int np = phase_usage_.num_phases;
        if (dkrds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i);
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalKr(s + np*i, kr + np*i);
            }
        }
    }




    /// Capillary pressure.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::capPress(const int n,
                                           const double* s,
                                           const int* cells,
                                           double* pc,
                                           double* dpcds) const
    {
        ASSERT (cells != 0);

        const int np = phase_usage_.num_phases;
        if (dpcds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i);
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).evalPc(s + np*i, pc + np*i);
            }
        }
    }




    /// Obtain the range of allowable saturation values.
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::satRange(const int n,
                                           const int* cells,
                                           double* smin,
                                           double* smax) const
    {
        ASSERT (cells != 0);

        const int np = phase_usage_.num_phases;
        for (int i = 0; i < n; ++i) {
            for (int p = 0; p < np; ++p) {
                smin[np*i + p] = funcForCell(cells[i]).smin_[p];
                smax[np*i + p] = funcForCell(cells[i]).smax_[p];
            }
        }
    }


    // Map the cell number to the correct function set.
    template <class SatFuncSet>
    const typename SaturationPropsFromDeck<SatFuncSet>::Funcs&
    SaturationPropsFromDeck<SatFuncSet>::funcForCell(const int cell) const
    {
        return cell_to_func_.empty() ? satfuncset_[0] : satfuncset_[cell_to_func_[cell]];
    }



} // namespace Opm

#endif // OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
