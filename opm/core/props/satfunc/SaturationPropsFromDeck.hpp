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

#include <opm/core/props/satfunc/SaturationPropsInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/satfunc/SatFuncMultiplexer.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <vector>

struct UnstructuredGrid;

namespace Opm
{
    /// Interface to saturation functions from deck.
    class SaturationPropsFromDeck : public SaturationPropsInterface
    {
    public:
        /// Default constructor.
        inline SaturationPropsFromDeck();

        /// Initialize from deck and grid.
        /// \param[in]  deck     Deck input parser
        /// \param[in]  grid     Grid to which property object applies, needed for the
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        inline void init(Opm::DeckConstPtr deck,
                         Opm::EclipseStateConstPtr eclipseState,
                         const UnstructuredGrid& grid);

        /// Initialize from deck and grid.
        /// \param[in]  deck     Deck input parser
        /// \param[in]  number_of_cells The number of cells of the grid to which property
        ///                             object applies, needed for the
        ///                             mapping from cell indices (typically from a processed
        ///                             grid) to logical cartesian indices consistent with the
        ///                             deck.
        /// \param[in]  global_cell     The mapping from local cell indices of the grid to
        ///                             global cell indices used in the deck.
        /// \param[in]  begin_cell_centroids Pointer to the first cell_centroid of the grid.
        /// \param[in]  dimensions      The dimensions of the grid. 
        template<class T>
        inline void init(Opm::DeckConstPtr deck,
                         Opm::EclipseStateConstPtr eclipseState,
                         int number_of_cells,
                         const int* global_cell,
                         const T& begin_cell_centroids,
                         int dimensions);

        /// \return   P, the number of phases.
        inline int numPhases() const;

        /// Relative permeability.
        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
        inline void relperm(const int n,
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
        inline void capPress(const int n,
                             const double* s,
                             const int* cells,
                             double* pc,
                             double* dpcds) const;

        /// Obtain the range of allowable saturation values.
        /// \param[in]  n      Number of data points.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        inline void satRange(const int n,
                             const int* cells,
                             double* smin,
                             double* smax) const;

        /// Update saturation state for the hysteresis tracking 
        /// \param[in]  n      Number of data points. 
        /// \param[in]  s      Array of nP saturation values.             
        inline void updateSatHyst(const int n,
                                  const int* cells,
                                  const double* s);

        /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
        /// \param[in]     cell  Cell index. 
        /// \param[in]     pcow  P_oil - P_water.
        /// \param[in/out] swat  Water saturation. / Possibly modified Water saturation.        
        inline void swatInitScaling(const int cell, 
                                    const double pcow, 
                                    double & swat);

    private:
        // internal helper method for satRange()
        template <class SaturationFunction>
        void satRange_(const SaturationFunction& satFunc,
                       const int cellIdx,
                       const int* cells,
                       double* smin,
                       double* smax) const;

        PhaseUsage phase_usage_;
        typedef Opm::SatFuncMultiplexer SatFunc;
        std::vector<SatFunc> satfunc_;
        std::vector<int> cell_to_func_; // = SATNUM - 1
        std::vector<int> cell_to_func_imb_;

        bool do_eps_;  // ENDSCALE is active
        bool do_3pt_;  // SCALECRS: YES~true  NO~false
        bool do_hyst_;  // Keywords ISWL etc detected     
        std::vector<EPSTransforms> eps_transf_;
        std::vector<EPSTransforms> eps_transf_hyst_;
        std::vector<SatHyst> sat_hyst_;

        template<class T>
        void initEPS(Opm::DeckConstPtr deck,
                     Opm::EclipseStateConstPtr eclipseState,
                     int number_of_cells,
                     const int* global_cell,
                     const T& begin_cell_centroids,
                     int dimensions,
                     const std::vector<std::string>& eps_kw,
                     std::vector<EPSTransforms>& eps_transf);
        template<class T>
        void initEPSKey(Opm::DeckConstPtr deck,
                        Opm::EclipseStateConstPtr eclipseState,
                        int number_of_cells,
                        const int* global_cell,
                        const T& begin_cell_centroids,
                        int dimensions,
                        const std::string& keyword,
                        std::vector<double>& scaleparam);
        void initEPSParam(const int cell, 
                          EPSTransforms::Transform& data,
                          const bool oil,
                          const double sl_tab,
                          const double scr_tab,
                          const double su_tab,
                          const double sxcr_tab,
                          const double s0_tab,
                          const double krsr_tab,
                          const double krmax_tab,
                          const double pcmax_tab,
                          const std::vector<double>& sl,
                          const std::vector<double>& scr,
                          const std::vector<double>& su,
                          const std::vector<double>& sxcr,
                          const std::vector<double>& s0,
                          const std::vector<double>& krsr,
                          const std::vector<double>& krmax,
                          const std::vector<double>& pcmax);

        bool columnIsMasked_(Opm::DeckConstPtr deck,
                             const std::string& keywordName,
                             int columnIdx)
        { return deck->getKeyword(keywordName)->getRecord(columnIdx)->getItem(0)->getSIDouble(0) != -1.0; }
    };



} // namespace Opm


#include <opm/core/props/satfunc/SaturationPropsFromDeck_impl.hpp>


#endif // OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
