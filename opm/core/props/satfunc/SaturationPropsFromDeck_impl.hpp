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

        // Saturation table scaling
        do_eps_ = false;
        do_3pt_ = false;
        if (deck.hasField("ENDSCALE")) {
            if (!phase_usage_.phase_used[Aqua] || !phase_usage_.phase_used[Liquid] || phase_usage_.phase_used[Vapour]) {
                THROW("Currently endpoint-scaling limited to oil-water systems without gas.");
            }
            if (deck.getENDSCALE().dir_switch_ != std::string("NODIR")) {
                THROW("SaturationPropsFromDeck::init()   --  ENDSCALE: Currently only 'NODIR' accepted.");
            }
            if (deck.getENDSCALE().revers_switch_ != std::string("REVERS")) {
                THROW("SaturationPropsFromDeck::init()   --  ENDSCALE: Currently only 'REVERS' accepted.");
            }
            if (deck.hasField("SCALECRS")) {
                if (deck.getSCALECRS().scalecrs_ == std::string("YES")) {
                    do_3pt_ = true;
                }
            }
            do_eps_ = true;
            initEPS(deck, grid, std::string("SWCR"), eps_.swcr_);
            initEPS(deck, grid, std::string("SWL"), eps_.swl_);
            initEPS(deck, grid, std::string("SWU"), eps_.swu_);
            initEPS(deck, grid, std::string("SOWCR"), eps_.sowcr_);
            initEPS(deck, grid, std::string("KRW"), eps_.krw_);
            initEPS(deck, grid, std::string("KRWR"), eps_.krwr_);
            initEPS(deck, grid, std::string("KRO"), eps_.kro_);
            initEPS(deck, grid, std::string("KRORW"), eps_.krorw_);
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
                if (do_eps_) {
                   relpermEPS(s + np*i, cells[i], kr + np*i, dkrds + np*np*i);
                } else {
                   funcForCell(cells[i]).evalKrDeriv(s + np*i, kr + np*i, dkrds + np*np*i);
                }
            }
        } else {
// #pragma omp parallel for
            for (int i = 0; i < n; ++i) {
                if (do_eps_) {
                   relpermEPS(s + np*i, cells[i], kr + np*i);
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

    // Initialize saturation scaling parameter
    template <class SatFuncSet>
    void SaturationPropsFromDeck<SatFuncSet>::initEPS(const EclipseGridParser& deck,
                                                      const UnstructuredGrid& grid,
                                                      const std::string& keyword,
                                                      std::vector<double>& scaleparam)
    {
        bool useKeyword = deck.hasField(keyword);  
        bool hasENPTVD = deck.hasField("ENPTVD");
        bool hasENKRVD = deck.hasField("ENKRVD");
        int itab = 0;
        std::vector<std::vector<std::vector<double> > > table_dummy;
        std::vector<std::vector<std::vector<double> > >& table = table_dummy;

        // Active keyword assigned default values for each cell (in case of possible box-wise assignment)
        int phase_pos_aqua = phase_usage_.phase_pos[BlackoilPhases::Aqua];
        if (keyword[0] == 'S' && (useKeyword || hasENPTVD)) {
            if (keyword == std::string("SWL")) {
                if (useKeyword || deck.getENPTVD().mask_[0]) {
                    itab = 1;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smin_[phase_pos_aqua];
                }
            } else if (keyword == std::string("SWCR")) {
                if (useKeyword || deck.getENPTVD().mask_[1]) {
                    itab = 2;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).swcr_;
                }
            } else if (keyword == std::string("SWU")) {
                if (useKeyword || deck.getENPTVD().mask_[2]) {
                    itab = 3;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).smax_[phase_pos_aqua];
                }
            } else if (keyword == std::string("SOWCR")) {
                if (useKeyword || deck.getENPTVD().mask_[3]) {
                    itab = 4;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).sowcr_;
                }
            }else {
                THROW(" -- unknown keyword: '" << keyword << "'");
            }
            if (!useKeyword && itab > 0) {
                table = deck.getENPTVD().table_;
            } 
        } else if (keyword[0] == 'K' && (useKeyword || hasENKRVD)) {    
            if (keyword == std::string("KRW")) {
                if (useKeyword || deck.getENKRVD().mask_[0]) {
                    itab = 1;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwmax_;
                }
            } else if (keyword == std::string("KRO")) {
                if (useKeyword || deck.getENKRVD().mask_[1]) {
                    itab = 2;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).kromax_;
                }
            } else if (keyword == std::string("KRWR")) {
                if (useKeyword || deck.getENKRVD().mask_[2]) {
                    itab = 3;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krwr_;
                }
            } else if (keyword == std::string("KRORW")) {
                if (useKeyword || deck.getENKRVD().mask_[3]) {
                    itab = 4;
                    scaleparam.resize(grid.number_of_cells);
                    for (int i=0; i<grid.number_of_cells; ++i)
                        scaleparam[i] = funcForCell(i).krorw_;
                }
            } else {
                THROW(" -- unknown keyword: '" << keyword << "'");
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
    void SaturationPropsFromDeck<SatFuncSet>::relpermEPS(const double *s, const int cell, double *kr, double *dkrds) const
    {
       const int wpos = phase_usage_.phase_pos[BlackoilPhases::Aqua];
       const int opos = phase_usage_.phase_pos[BlackoilPhases::Liquid];
       double ss[PhaseUsage::MaxNumPhases];
       
       if (do_3pt_) { // Three-point scaling
           // Transforms for water saturation 
           if (eps_.swcr_.empty() && eps_.swu_.empty()) {
               ss[wpos] = s[wpos]; 
           } else {
               double s_r = 1.0-funcForCell(cell).sowcr_;
               double sr = eps_.sowcr_.empty() ? s_r : 1.0-eps_.sowcr_[cell];
               if (s[wpos] <= sr) {
                   double sw_cr = funcForCell(cell).swcr_;
                   double swcr = eps_.swcr_.empty() ? sw_cr : eps_.swcr_[cell];
                   ss[wpos] = (s[wpos] <= swcr) ? sw_cr : sw_cr+(s[wpos]-swcr)*(s_r-sw_cr)/(sr-swcr);
               } else {
                   double sw_max = funcForCell(cell).smax_[wpos];
                   double swmax = eps_.swu_.empty() ? sw_max : eps_.swu_[cell];
                   ss[wpos] = (s[wpos] >= swmax) ? sw_max : s_r+(s[wpos]-sr)*(sw_max-s_r)/(swmax-sr);
               }
           }
           // Transforms for oil saturation 
           if (eps_.sowcr_.empty() && eps_.swl_.empty()) {
               ss[opos] = s[opos]; 
           } else {
               double s_r = 1.0-funcForCell(cell).swcr_;
               double sr = eps_.swcr_.empty() ? s_r : 1.0-eps_.swcr_[cell];
               if (s[opos] <= sr) {
                   double sow_cr = funcForCell(cell).sowcr_;
                   double sowcr = eps_.sowcr_.empty() ? sow_cr : eps_.sowcr_[cell];
                   ss[opos] = (s[opos] <= sowcr) ? sow_cr : sow_cr+(s[opos]-sowcr)*(s_r-sow_cr)/(sr-sowcr);
               } else {
                   double sow_max = funcForCell(cell).smax_[opos];
                   double sowmax = eps_.swl_.empty() ? sow_max : (1.0-eps_.swl_[cell]);
                   ss[opos] = (s[opos] >= sowmax) ? sow_max : s_r+(s[opos]-sr)*(sow_max-s_r)/(sowmax-sr);
               }
           }
       } else { // Two-point scaling
           // Transforms for water saturation 
           if (eps_.swcr_.empty() && eps_.swu_.empty()) {
               ss[wpos] = s[wpos]; 
           } else {
               double sw_cr = funcForCell(cell).swcr_;
               double swcr = eps_.swcr_.empty() ? sw_cr : eps_.swcr_[cell];
               if (s[wpos] <= swcr) {
                   ss[wpos] = sw_cr;
               } else {
                   double sw_max = funcForCell(cell).smax_[wpos];
                   double swmax = eps_.swu_.empty() ? sw_max : eps_.swu_[cell];
                   ss[wpos] = (s[wpos] >= swmax) ? sw_max : sw_cr + (s[wpos]-swcr)*(sw_max-sw_cr)/(swmax-swcr);
               }
           }
           // Transforms for oil saturation 
           if (eps_.sowcr_.empty() && eps_.swl_.empty()) {
               ss[opos] = s[opos]; 
           } else {
               double sow_cr = funcForCell(cell).sowcr_;
               double socr = eps_.sowcr_.empty() ? sow_cr : eps_.sowcr_[cell];
               if (s[opos] <= socr) {
                   ss[opos] = sow_cr;
               } else {
                   double sow_max = funcForCell(cell).smax_[opos];
                   double sowmax = eps_.swl_.empty() ? sow_max : (1.0-eps_.swl_[cell]);
                   ss[opos] = (s[opos] >= sowmax) ? sow_max : sow_cr + (s[opos]-socr) *(sow_max-sow_cr)/(sowmax-socr);
               }
           }
       }

       // Evaluation of relperms
       if (dkrds) {
           THROW("Relperm derivatives not yet available in combination with end point scaling ...");
           funcForCell(cell).evalKrDeriv(ss, kr, dkrds);
       } else {
           // Assume: sw_cr -> krw=0     sw_max -> krw=<max water relperm>
           //         sow_cr -> kro=0    sow_max -> kro=<max oil relperm>
           funcForCell(cell).evalKr(ss, kr);
       }  

       // Scaling of relperms values
       //  - Water
       if (eps_.krw_.empty() && eps_.krwr_.empty()) { // No value scaling
       } else if (eps_.krwr_.empty()) { // Two-point
           kr[wpos] *= (eps_.krw_[cell]/funcForCell(cell).krwmax_);
       } else {
           double swcr = eps_.swcr_.empty() ? funcForCell(cell).swcr_ : eps_.swcr_[cell];
           double swmax = eps_.swu_.empty() ? funcForCell(cell).smax_[wpos] : eps_.swu_[cell];
           double sr;
           if (do_3pt_) {
               sr = eps_.sowcr_.empty() ? 1.0-funcForCell(cell).sowcr_ : 1.0-eps_.sowcr_[cell];
           } else {
               double sw_cr = funcForCell(cell).swcr_;
               double sw_max = funcForCell(cell).smax_[wpos];
               double s_r = 1.0-funcForCell(cell).sowcr_;
               sr = swcr + (s_r-sw_cr)*(swmax-swcr)/(sw_max-sw_cr);
           }          
           if (s[wpos] <= swcr) {
               kr[wpos] = 0.0;
           } else if (sr > swmax-1.0e-6) {
               if (do_3pt_) { //Ignore krw and do two-point?
                   kr[wpos] *= eps_.krwr_[cell]/funcForCell(cell).krwr_;
               } else if (!eps_.kro_.empty()){ //Ignore krwr and do two-point
                   kr[wpos] *= eps_.krw_[cell]/funcForCell(cell).krwmax_;
               }
           } else if (s[wpos] <= sr) {
               kr[wpos] *= eps_.krwr_[cell]/funcForCell(cell).krwr_;
           } else if (s[wpos] <= swmax) {
               double krw_max = funcForCell(cell).krwmax_;
               double krw = eps_.krw_.empty() ? krw_max : eps_.krw_[cell];
               double krw_r = funcForCell(cell).krwr_;
               double krwr = eps_.krwr_.empty() ? krw_r : eps_.krwr_[cell];
               if (std::fabs(krw_max- krw_r) > 1.0e-6) {
                   kr[wpos] = krwr + (kr[wpos]-krw_r)*(krw-krwr)/(krw_max-krw_r);
               } else {
                   kr[wpos] = krwr + (krw-krwr)*(s[wpos]-sr)/(swmax-sr);
               }
           } else {
               kr[wpos] = eps_.krw_.empty() ? funcForCell(cell).krwmax_ : eps_.krw_[cell];
           }
       }
       
       //  - Oil
       if (eps_.kro_.empty() && eps_.krorw_.empty()) { // No value scaling
       } else if (eps_.krorw_.empty()) { // Two-point scaling
           kr[opos] *= (eps_.kro_[cell]/funcForCell(cell).kromax_);
       } else {
           double sowcr = eps_.sowcr_.empty() ? funcForCell(cell).sowcr_ : eps_.sowcr_[cell];
           double sowmax = eps_.swl_.empty() ? funcForCell(cell).smax_[opos] : 1.0-eps_.swl_[cell];
           double sr;
           if (do_3pt_) {
               sr = eps_.swcr_.empty() ? 1.0-funcForCell(cell).swcr_ : 1.0-eps_.swcr_[cell];
           } else {
               double sow_cr = funcForCell(cell).sowcr_;
               double sow_max = funcForCell(cell).smax_[opos];
               double s_r = 1.0-funcForCell(cell).swcr_;
               sr = sowcr + (s_r-sow_cr)*(sowmax-sowcr)/(sow_max-sow_cr);
           }
           if (s[opos] <= sowcr) {
               kr[opos] = 0.0;
           } else if (sr > sowmax-1.0e-6) {
               if (do_3pt_) { //Ignore kro and do two-point?
                   kr[opos] *= eps_.krorw_[cell]/funcForCell(cell).krorw_;
               } else if (!eps_.kro_.empty()){ //Ignore krowr and do two-point
                   kr[opos] *= eps_.kro_[cell]/funcForCell(cell).kromax_;
               }
           } else if (s[opos] <= sr) {
               kr[opos] *= eps_.krorw_[cell]/funcForCell(cell).krorw_;
           } else if (s[opos] <= sowmax) {
               double kro_max = funcForCell(cell).kromax_;
               double kro = eps_.kro_.empty() ? kro_max : eps_.kro_[cell];
               double kro_rw = funcForCell(cell).krorw_;
               double krorw = eps_.krorw_[cell];
               if (std::fabs(kro_max- kro_rw) > 1.0e-6) {
                   kr[opos] = krorw + (kr[opos]- kro_rw)*(kro-krorw)/(kro_max- kro_rw);
               } else {
                   kr[opos] = krorw + (kro-krorw)*(s[opos]-sr)/(sowmax-sr);
               }
           } else {
               kr[opos] = eps_.kro_.empty() ? funcForCell(cell).kromax_ : eps_.kro_[cell];
           }
       }
    }

} // namespace Opm

#endif // OPM_SATURATIONPROPSFROMDECK_IMPL_HEADER_INCLUDED
