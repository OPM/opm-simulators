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
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/grid.h>

#include <iostream>

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
            OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  oil phase must be active.");
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
                OPM_THROW(std::runtime_error, "Found " << num_tables << " SWOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
        }
        if (phase_usage_.phase_used[Vapour]) {
            const SGOF::table_t& sgof_table = deck.getSGOF().sgof_;
            int num_sgof_tables = sgof_table.size();
            if (num_sgof_tables < satfuncs_expected) {
                OPM_THROW(std::runtime_error, "Found " << num_tables << " SGOF tables, SATNUM specifies at least " << satfuncs_expected);
            }
            if (num_tables == Uninitialized) {
                num_tables = num_sgof_tables;
            } else if (num_tables != num_sgof_tables) {
                OPM_THROW(std::runtime_error, "Inconsistent number of tables in SWOF and SGOF.");
            }
        }

        // Initialize tables.
        satfuncset_.resize(num_tables);
        for (int table = 0; table < num_tables; ++table) {
            satfuncset_[table].init(deck, table, phase_usage_, samples);
        }

        // Saturation table scaling
        do_eps_ = false;
        do_3pt_ = false;
        if (deck.hasField("ENDSCALE")) {
            //if (!phase_usage_.phase_used[Aqua] || !phase_usage_.phase_used[Liquid] || phase_usage_.phase_used[Vapour]) {
            //    OPM_THROW(std::runtime_error, "Currently endpoint-scaling limited to oil-water systems without gas.");
            //}
            if (deck.getENDSCALE().dir_switch_ != std::string("NODIR")) {
                OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  ENDSCALE: Currently only 'NODIR' accepted.");
            }
            if (deck.getENDSCALE().revers_switch_ != std::string("REVERS")) {
                OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  ENDSCALE: Currently only 'REVERS' accepted.");
            }
            if (deck.hasField("SCALECRS")) {
                if (deck.getSCALECRS().scalecrs_ == std::string("YES")) {
                    do_3pt_ = true;
                }
            }
            do_eps_ = true;
            
            initEPS(deck, grid);
            
            // For now, a primitive detection of hysteresis. TODO: SATOPTS HYSTER/ and EHYSTR
            do_hyst_ = deck.hasField("ISWL") || deck.hasField("ISWU") || deck.hasField("ISWCR") || deck.hasField("ISGL") ||
                       deck.hasField("ISGU") || deck.hasField("ISGCR") ||  deck.hasField("ISOWCR") || deck.hasField("ISOGCR");
            if (do_hyst_) {
                if (deck.hasField("KRW") || deck.hasField("KRG") || deck.hasField("KRO") || deck.hasField("KRWR") || 
                    deck.hasField("KRGR") || deck.hasField("KRORW") || deck.hasField("KRORG") ||
                    deck.hasField("IKRW") || deck.hasField("IKRG") || deck.hasField("IKRO") || deck.hasField("IKRWR") || 
                    deck.hasField("IKRGR") || deck.hasField("IKRORW") || deck.hasField("IKRORG") ) {
                    OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  ENDSCALE: Currently hysteresis and relperm value scaling can not be combined.");
                }
                initEPSHyst(deck, grid);
            }

            //OPM_THROW(std::runtime_error, "SaturationPropsFromDeck::init()   --  ENDSCALE: Under construction ...");
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
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (dkrds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_hyst_) {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
                } else if (do_eps_) {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i);
                }
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_hyst_) {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
                } else if (do_eps_) {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalKr(s + np*i, kr + np*i);
                }
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
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (dpcds) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_eps_) {
                   funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalPcDeriv(s + np*i, pc + np*i, dpcds + np*np*i);
                }
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {         
                if (do_eps_) {
                   funcForCell(cells[i]).evalPc(s + np*i, pc + np*i, &(eps_transf_[cells[i]]));
                } else {
                   funcForCell(cells[i]).evalPc(s + np*i, pc + np*i);
                }
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
        assert(cells != 0);
        const int np = phase_usage_.num_phases;
       
        if (do_eps_) {
            const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
            const int opos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
            const int gpos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
            for (int i = 0; i < n; ++i) {
                smin[np*i + opos] = 1.0;
                smax[np*i + opos] = 1.0;
                if (phase_usage_.phase_used[Aqua]) {
                    smin[np*i + wpos] = eps_transf_[cells[i]].wat.doNotScale ? eps_transf_[cells[i]].wat.smin: funcForCell(cells[i]).smin_[wpos];
                    smax[np*i + wpos] = eps_transf_[cells[i]].wat.doNotScale ? eps_transf_[cells[i]].wat.smax: funcForCell(cells[i]).smax_[wpos];
                    smin[np*i + opos] -= smax[np*i + wpos];
                    smax[np*i + opos] -= smin[np*i + wpos];
                }  
                if (phase_usage_.phase_used[Vapour]) {
                    smin[np*i + gpos] = eps_transf_[cells[i]].wat.doNotScale ? eps_transf_[cells[i]].gas.smin: funcForCell(cells[i]).smin_[gpos];
                    smax[np*i + gpos] = eps_transf_[cells[i]].wat.doNotScale ? eps_transf_[cells[i]].gas.smax: funcForCell(cells[i]).smax_[gpos];
                    smin[np*i + opos] -= smax[np*i + gpos];
                    smax[np*i + opos] -= smin[np*i + gpos];
                }
                if (phase_usage_.phase_used[Vapour] && phase_usage_.phase_used[Aqua]) {
                    smin[np*i + opos] = std::max(0.0,smin[np*i + opos]);
                }
            }
        } else {
            for (int i = 0; i < n; ++i) {
                for (int p = 0; p < np; ++p) {
                    smin[np*i + p] = funcForCell(cells[i]).smin_[p];
                    smax[np*i + p] = funcForCell(cells[i]).smax_[p];
                }
            }
        }
    }

        
    /// Update saturation state for the hysteresis tracking 
    /// \param[in]  n      Number of data points. 
    /// \param[in]  s      Array of nP saturation values.
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::updateSatHyst(const int n,
                                                            const int* cells,
                                                            const double* s)
    {        
        assert(cells != 0);

        const int np = phase_usage_.num_phases;
        if (do_hyst_) {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                funcForCell(cells[i]).updateSatHyst(s + np*i, &(eps_transf_[cells[i]]), &(eps_transf_hyst_[cells[i]]), &(sat_hyst_[cells[i]]));
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

    // Initialize saturation scaling parameters
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPS(const EclipseGridParser& deck,
                                                      const UnstructuredGrid& grid)
    {
        std::vector<double> swl, swcr, swu, sgl, sgcr, sgu, sowcr, sogcr;
        std::vector<double> krw, krg, kro, krwr, krgr, krorw, krorg;
        // Initialize saturation scaling parameter
        initEPSKey(deck, grid, std::string("SWL"),   swl);
        initEPSKey(deck, grid, std::string("SWU"),   swu);
        initEPSKey(deck, grid, std::string("SWCR"),  swcr);
        initEPSKey(deck, grid, std::string("SGL"),   sgl);
        initEPSKey(deck, grid, std::string("SGU"),   sgu);
        initEPSKey(deck, grid, std::string("SGCR"),  sgcr);
        initEPSKey(deck, grid, std::string("SOWCR"), sowcr);
        initEPSKey(deck, grid, std::string("SOGCR"), sogcr);
        initEPSKey(deck, grid, std::string("KRW"),   krw);
        initEPSKey(deck, grid, std::string("KRG"),   krg);
        initEPSKey(deck, grid, std::string("KRO"),   kro);
        initEPSKey(deck, grid, std::string("KRWR"),  krwr);
        initEPSKey(deck, grid, std::string("KRGR"),  krgr);
        initEPSKey(deck, grid, std::string("KRORW"), krorw);
        initEPSKey(deck, grid, std::string("KRORG"), krorg);

        eps_transf_.resize(grid.number_of_cells);

        const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        const int gpos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        const bool oilWater = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && !phase_usage_.phase_used[Vapour];
        const bool oilGas = !phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];
        const bool threephase = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];

        for (int cell = 0; cell < grid.number_of_cells; ++cell) {
            if (oilWater) {
                // ### krw
                initEPSParam(cell, eps_transf_[cell].wat, false, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_, funcForCell(cell).smax_[wpos],
                  funcForCell(cell).sowcr_, -1.0, funcForCell(cell).krwr_, funcForCell(cell).krwmax_, swl, swcr, swu, sowcr, sgl, krwr, krw);
                // ### krow
                initEPSParam(cell, eps_transf_[cell].watoil, true, 0.0, funcForCell(cell).sowcr_, funcForCell(cell).smin_[wpos],
                  funcForCell(cell).swcr_, -1.0, funcForCell(cell).krorw_, funcForCell(cell).kromax_, swl, sowcr, swl, swcr, sgl, krorw, kro);
            } else if (oilGas) {
                // ### krg
                initEPSParam(cell, eps_transf_[cell].gas, false, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_, funcForCell(cell).smax_[gpos],
                  funcForCell(cell).sogcr_, -1.0, funcForCell(cell).krgr_, funcForCell(cell).krgmax_, sgl, sgcr, sgu, sogcr, swl, krgr, krg);
                // ### krog
                initEPSParam(cell, eps_transf_[cell].gasoil, true, 0.0, funcForCell(cell).sogcr_, funcForCell(cell).smin_[gpos],
                  funcForCell(cell).sgcr_, -1.0, funcForCell(cell).krorg_, funcForCell(cell).kromax_, sgl, sogcr, sgl, sgcr, swl, krorg, kro);
            } else if (threephase) {
                // ### krw
                initEPSParam(cell, eps_transf_[cell].wat, false, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_, funcForCell(cell).smax_[wpos], funcForCell(cell).sowcr_,
                  funcForCell(cell).smin_[gpos], funcForCell(cell).krwr_, funcForCell(cell).krwmax_, swl, swcr, swu, sowcr, sgl, krwr, krw);
                // ### krow
                initEPSParam(cell, eps_transf_[cell].watoil, true, 0.0, funcForCell(cell).sowcr_, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_,
                  funcForCell(cell).smin_[gpos], funcForCell(cell).krorw_, funcForCell(cell).kromax_, swl, sowcr, swl, swcr, sgl, krorw, kro);
                // ### krg
                initEPSParam(cell, eps_transf_[cell].gas, false, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_, funcForCell(cell).smax_[gpos], funcForCell(cell).sogcr_,
                  funcForCell(cell).smin_[wpos], funcForCell(cell).krgr_, funcForCell(cell).krgmax_, sgl, sgcr, sgu, sogcr, swl, krgr, krg);
                // ### krog
                initEPSParam(cell, eps_transf_[cell].gasoil, true, 0.0, funcForCell(cell).sogcr_, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_,
                  funcForCell(cell).smin_[wpos], funcForCell(cell).krorg_, funcForCell(cell).kromax_, sgl, sogcr, sgl, sgcr, swl, krorg, kro);
            }
        }
    }


    // Initialize hysteresis saturation scaling parameters
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPSHyst(const EclipseGridParser& deck,
                                                          const UnstructuredGrid& grid)
    {
        std::vector<double> iswl, iswcr, iswu, isgl, isgcr, isgu, isowcr, isogcr;
        std::vector<double> ikrw, ikrg, ikro, ikrwr, ikrgr, ikrorw, ikrorg;
        // Initialize hysteresis saturation scaling parameters
        initEPSKey(deck, grid, std::string("ISWL"),   iswl);
        initEPSKey(deck, grid, std::string("ISWU"),   iswu);
        initEPSKey(deck, grid, std::string("ISWCR"),  iswcr);
        initEPSKey(deck, grid, std::string("ISGL"),   isgl);
        initEPSKey(deck, grid, std::string("ISGU"),   isgu);
        initEPSKey(deck, grid, std::string("ISGCR"),  isgcr);
        initEPSKey(deck, grid, std::string("ISOWCR"), isowcr);
        initEPSKey(deck, grid, std::string("ISOGCR"), isogcr);
        initEPSKey(deck, grid, std::string("IKRW"),   ikrw);
        initEPSKey(deck, grid, std::string("IKRG"),   ikrg);
        initEPSKey(deck, grid, std::string("IKRO"),   ikro);
        initEPSKey(deck, grid, std::string("IKRWR"),  ikrwr);
        initEPSKey(deck, grid, std::string("IKRGR"),  ikrgr);
        initEPSKey(deck, grid, std::string("IKRORW"), ikrorw);
        initEPSKey(deck, grid, std::string("IKRORG"), ikrorg);
        

        eps_transf_hyst_.resize(grid.number_of_cells);
        sat_hyst_.resize(grid.number_of_cells);

        const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        const int gpos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        const bool oilWater = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && !phase_usage_.phase_used[Vapour];
        const bool oilGas = !phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];
        const bool threephase = phase_usage_.phase_used[Aqua] && phase_usage_.phase_used[Liquid] && phase_usage_.phase_used[Vapour];

        for (int cell = 0; cell < grid.number_of_cells; ++cell) {
            if (oilWater) {
                // ### krw
                initEPSParam(cell, eps_transf_hyst_[cell].wat, false, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_, funcForCell(cell).smax_[wpos],
                  funcForCell(cell).sowcr_, -1.0, funcForCell(cell).krwr_, funcForCell(cell).krwmax_, iswl, iswcr, iswu, isowcr, isgl, ikrwr, ikrw);
                // ### krow
                initEPSParam(cell, eps_transf_hyst_[cell].watoil, true, 0.0, funcForCell(cell).sowcr_, funcForCell(cell).smin_[wpos],
                  funcForCell(cell).swcr_, -1.0, funcForCell(cell).krorw_, funcForCell(cell).kromax_, iswl, isowcr, iswl, iswcr, isgl, ikrorw, ikro);
            } else if (oilGas) {
                // ### krg
                initEPSParam(cell, eps_transf_hyst_[cell].gas, false, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_, funcForCell(cell).smax_[gpos],
                  funcForCell(cell).sogcr_, -1.0, funcForCell(cell).krgr_, funcForCell(cell).krgmax_, isgl, isgcr, isgu, isogcr, iswl, ikrgr, ikrg);
                // ### krog
                initEPSParam(cell, eps_transf_hyst_[cell].gasoil, true, 0.0, funcForCell(cell).sogcr_, funcForCell(cell).smin_[gpos],
                  funcForCell(cell).sgcr_, -1.0, funcForCell(cell).krorg_, funcForCell(cell).kromax_, isgl, isogcr, isgl, isgcr, iswl, ikrorg, ikro);
            } else if (threephase) {
                // ### krw
                initEPSParam(cell, eps_transf_hyst_[cell].wat, false, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_, funcForCell(cell).smax_[wpos], funcForCell(cell).sowcr_,
                  funcForCell(cell).smin_[gpos], funcForCell(cell).krwr_, funcForCell(cell).krwmax_, iswl, iswcr, iswu, isowcr, isgl, ikrwr, ikrw);
                // ### krow
                initEPSParam(cell, eps_transf_hyst_[cell].watoil, true, 0.0, funcForCell(cell).sowcr_, funcForCell(cell).smin_[wpos], funcForCell(cell).swcr_,
                  funcForCell(cell).smin_[gpos], funcForCell(cell).krorw_, funcForCell(cell).kromax_, iswl, isowcr, iswl, iswcr, isgl, ikrorw, ikro);
                // ### krg
                initEPSParam(cell, eps_transf_hyst_[cell].gas, false, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_, funcForCell(cell).smax_[gpos], funcForCell(cell).sogcr_,
                  funcForCell(cell).smin_[wpos], funcForCell(cell).krgr_, funcForCell(cell).krgmax_, isgl, isgcr, isgu, isogcr, iswl, ikrgr, ikrg);
                // ### krog
                initEPSParam(cell, eps_transf_hyst_[cell].gasoil, true, 0.0, funcForCell(cell).sogcr_, funcForCell(cell).smin_[gpos], funcForCell(cell).sgcr_,
                  funcForCell(cell).smin_[wpos], funcForCell(cell).krorg_, funcForCell(cell).kromax_, isgl, isogcr, isgl, isgcr, iswl, ikrorg, ikro);
            }
        }
    }

    // Initialize saturation scaling parameter
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPSKey(const EclipseGridParser& deck,
                                                         const UnstructuredGrid& grid,
                                                         const std::string& keyword,
                                                         std::vector<double>& scaleparam)
    {
        const bool useAqua = phase_usage_.phase_used[Aqua];
        const bool useLiquid = phase_usage_.phase_used[Liquid];
        const bool useVapour = phase_usage_.phase_used[Vapour];
        bool useKeyword = deck.hasField(keyword);
        bool hasENPTVD = deck.hasField("ENPTVD");
        bool hasENKRVD = deck.hasField("ENKRVD");
        int itab = 0;
        std::vector<std::vector<std::vector<double> > > table_dummy;
        std::vector<std::vector<std::vector<double> > >& table = table_dummy;

        // Active keyword assigned default values for each cell (in case of possible box-wise assignment)
        int phase_pos_aqua = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        int phase_pos_vapour = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        if ((keyword[0] == 'S' && (useKeyword || hasENPTVD)) || (keyword[1] == 'S' && useKeyword) ) {
            if (keyword == std::string("SWL") || keyword == std::string("ISWL") ) {
                if (useAqua && (useKeyword || deck.getENPTVD().mask_[0])) {
                    itab = 1;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smin_[phase_pos_aqua];
                }
            } else if (keyword == std::string("SWCR") || keyword == std::string("ISWCR") ) {
                if (useAqua && (useKeyword || deck.getENPTVD().mask_[1])) {
                    itab = 2;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).swcr_;
                }
            } else if (keyword == std::string("SWU") || keyword == std::string("ISWU") ) {
                if (useAqua && (useKeyword || deck.getENPTVD().mask_[2])) {
                    itab = 3;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smax_[phase_pos_aqua];
                }
            } else if (keyword == std::string("SGL") || keyword == std::string("ISGL") ) {
                if (useVapour && (useKeyword || deck.getENPTVD().mask_[3])) {
                    itab = 4;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smin_[phase_pos_vapour];
                }
            } else if (keyword == std::string("SGCR") || keyword == std::string("ISGCR") ) {
                if (useVapour && (useKeyword || deck.getENPTVD().mask_[4])) {
                    itab = 5;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).sgcr_;
                }
            } else if (keyword == std::string("SGU") || keyword == std::string("ISGU") ) {
                if (useVapour && (useKeyword || deck.getENPTVD().mask_[5])) {
                    itab = 6;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smax_[phase_pos_vapour];
                }
            } else if (keyword == std::string("SOWCR") || keyword == std::string("ISOWCR") ) {
                if (useAqua && (useKeyword || deck.getENPTVD().mask_[6])) {
                    itab = 7;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).sowcr_;
                }
            } else if (keyword == std::string("SOGCR") || keyword == std::string("ISOGCR") ) {
                if (useVapour && (useKeyword || deck.getENPTVD().mask_[7])) {
                    itab = 8;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).sogcr_;
                }
            } else {
                OPM_THROW(std::runtime_error, " -- unknown keyword: '" << keyword << "'");
            }
            if (!useKeyword && itab > 0) {
                table = deck.getENPTVD().table_;
            }
        } else if ((keyword[0] == 'K' && (useKeyword || hasENKRVD)) || (keyword[1] == 'K' && useKeyword) ) {
            if (keyword == std::string("KRW") || keyword == std::string("IKRW") ) {
                if (useAqua && (useKeyword || deck.getENKRVD().mask_[0])) {
                    itab = 1;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwmax_;
                }
            } else if (keyword == std::string("KRG") || keyword == std::string("IKRG") ) {
                if (useVapour && (useKeyword || deck.getENKRVD().mask_[1])) {
                    itab = 2;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krgmax_;
                }
            } else if (keyword == std::string("KRO") || keyword == std::string("IKRO") ) {
                if (useLiquid && (useKeyword || deck.getENKRVD().mask_[2])) {
                    itab = 3;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).kromax_;
                }
            } else if (keyword == std::string("KRWR") || keyword == std::string("IKRWR") ) {
                if (useAqua && (useKeyword || deck.getENKRVD().mask_[3])) {
                    itab = 4;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwr_;
                }
            } else if (keyword == std::string("KRGR") || keyword == std::string("IKRGR") ) {
                if (useVapour && (useKeyword || deck.getENKRVD().mask_[4])) {
                    itab = 5;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krgr_;
                }
            } else if (keyword == std::string("KRORW") || keyword == std::string("IKRORW") ) {
                if (useAqua && (useKeyword || deck.getENKRVD().mask_[5])) {
                    itab = 6;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krorw_;
                }
            } else if (keyword == std::string("KRORG") || keyword == std::string("IKRORG") ) {
                if (useVapour && (useKeyword || deck.getENKRVD().mask_[6])) {
                    itab = 7;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krorg_;
                }
            } else {
                OPM_THROW(std::runtime_error, " -- unknown keyword: '" << keyword << "'");
            }
            if (!useKeyword && itab > 0) {
                table = deck.getENKRVD().table_;
            }
        }

        if (scaleparam.empty()) {
            return;
        } else if (useKeyword) {
            // Keyword values from deck
            std::cout << "--- Scaling parameter '" << keyword << "' assigned." << std::endl;
            const int* gc = grid.global_cell;
            const std::vector<double>& val = deck.getFloatingPointValue(keyword);
            for (int c = 0; c < int(scaleparam.size()); ++c) {
                const int deck_pos = (gc == NULL) ? c : gc[c];
                scaleparam[c] = val[deck_pos];
            }
        } else {
            std::cout << "--- Scaling parameter '" << keyword << "' assigned via ";
            if (keyword[0] == 'S')
                deck.getENPTVD().write(std::cout);
            else
                deck.getENKRVD().write(std::cout);
            const double* cc = grid.cell_centroids;
            const int dim = grid.dimensions;
            for (int cell = 0; cell < grid.number_of_cells; ++cell) {
                int jtab = cell_to_func_.empty() ? 0 : cell_to_func_[cell];
                if (table[itab][jtab][0] != -1.0) {
                    std::vector<double>& depth = table[0][jtab];
                    std::vector<double>& val = table[itab][jtab];
                    double zc = cc[dim*cell+dim-1];
                    if (zc >= depth.front() && zc <= depth.back()) { //don't want extrap outside depth interval
                        scaleparam[cell] = linearInterpolation(depth, val, zc);
                    }
                }
            }
        }
    }

    // Saturation scaling
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPSParam(const int cell,
                                                           EPSTransforms::Transform& data,
                                                           const bool oil,          // flag indicating krow/krog calculations
                                                           const double sl_tab,     // minimum saturation (for krow/krog calculations this is normally zero)
                                                           const double scr_tab,    // critical saturation
                                                           const double su_tab,     // maximum saturation (for krow/krog calculations this is minimum water/gas saturation)
                                                           const double sxcr_tab,   // second critical saturation (not used for 2pt scaling)
                                                           const double s0_tab,     // threephase complementary minimum saturation (-1.0 indicates 2-phase)
                                                           const double krsr_tab,   // relperm at displacing critical saturation
                                                           const double krmax_tab,  // relperm at maximum saturation
                                                           const std::vector<double>& sl,  // For krow/krog calculations this is not used
                                                           const std::vector<double>& scr,
                                                           const std::vector<double>& su,  // For krow/krog calculations this is SWL/SGL
                                                           const std::vector<double>& sxcr,
                                                           const std::vector<double>& s0,
                                                           const std::vector<double>& krsr,
                                                           const std::vector<double>& krmax)
    {
        if (scr.empty() && su.empty() && (sxcr.empty() || !do_3pt_) && s0.empty()) {
            data.doNotScale = true;
        } else {
            data.doNotScale = false;
            data.do_3pt = do_3pt_;
            double s_r;
            if (s0_tab < 0.0) { // 2phase
                s_r = 1.0-sxcr_tab;
                if (do_3pt_) data.sr = sxcr.empty() ? s_r : 1.0-sxcr[cell];
            } else { // 3phase
                s_r = 1.0-sxcr_tab-s0_tab;
                if (do_3pt_)data.sr = 1.0 - (sxcr.empty() ? sxcr_tab : sxcr[cell])
                                          - (s0.empty() ? s0_tab : s0[cell]);
            }
            data.scr = scr.empty() ? scr_tab : scr[cell];
            double s_max = su_tab;
            if (oil) {
                data.smin = sl_tab;
                if (s0_tab < 0.0) { // 2phase
                    s_max = 1.0 - su_tab;
                    data.smax = 1.0 - (su.empty() ? su_tab : su[cell]);
                } else { // 3phase
                    s_max = 1.0 - su_tab - s0_tab;
                    data.smax = 1.0 - (su.empty() ? su_tab : su[cell])
                                    - (s0.empty() ? s0_tab : s0[cell]);
                }
            } else {
                data.smin = sl.empty() ? sl_tab : sl[cell];
                data.smax = su.empty() ? su_tab : su[cell];
            }
            if (do_3pt_) {
                data.slope1 = (s_r-scr_tab)/(data.sr-data.scr);
                data.slope2 = (s_max-s_r)/(data.smax-data.sr);
            } else {
                data.slope2 = data.slope1 = (s_max-scr_tab)/(data.smax-data.scr);
                // Inv transform of tabulated critical displacing saturation to prepare for possible value scaling (krwr etc)
                data.sr = data.scr + (s_r-scr_tab)*(data.smax-data.scr)/(s_max-scr_tab);
            }
        }
        
        data.doKrMax = !krmax.empty();
        data.doKrCrit = !krsr.empty();
        data.doSatInterp = false;
        data.krsr = krsr.empty() ? krsr_tab : krsr[cell];
        data.krmax = krmax.empty() ? krmax_tab : krmax[cell];
        data.krSlopeCrit = data.krsr/krsr_tab;
        data.krSlopeMax = data.krmax/krmax_tab;
        if (data.doKrCrit) {
            if (data.sr > data.smax-1.0e-6) {
                //Ignore krsr and do two-point (one might consider combining krsr and krmax linearly between scr and smax ... )
                data.doKrCrit = false;
            } else if (std::fabs(krmax_tab- krsr_tab) > 1.0e-6) { // interpolate kr
                data.krSlopeMax = (data.krmax-data.krsr)/(krmax_tab-krsr_tab);
            } else { // interpolate sat
                data.doSatInterp = true;
                data.krSlopeMax = (data.krmax-data.krsr)/(data.smax-data.sr);
            }
        }

    }

} // namespace Opm

#endif // OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
