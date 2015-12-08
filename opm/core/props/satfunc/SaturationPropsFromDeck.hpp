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
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/grid.h>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    // Forward declaring the EclMaterialLawManager template.
    template <class ScalarT, int wettingPhaseIdxV, int nonWettingasPhaseIdxV, int gasPhaseIdxV>
    class ThreePhaseMaterialTraits;
    template <class Traits>
    class EclMaterialLawManager;


    /// Interface to saturation functions from deck.
    class SaturationPropsFromDeck : public SaturationPropsInterface
    {
    public:
        typedef Opm::ThreePhaseMaterialTraits<double,
                                              /*wettingPhaseIdx=*/BlackoilPhases::Aqua,
                                              /*nonWettingPhaseIdx=*/BlackoilPhases::Liquid,
                                              /*gasPhaseIdx=*/BlackoilPhases::Vapour> MaterialTraits;
        typedef Opm::EclMaterialLawManager<MaterialTraits> MaterialLawManager;

        /// Default constructor.
        SaturationPropsFromDeck();

        /// Initialize from a MaterialLawManager object.
        /// \param[in]  phaseUsage          Phase configuration
        /// \param[in]  materialLawManager  An initialized MaterialLawManager object
        void init(const PhaseUsage& phaseUsage,
                  std::shared_ptr<MaterialLawManager> materialLawManager);


        /// Initialize from deck and MaterialLawManager.
        /// \param[in]  deck                Input deck
        /// \param[in]  materialLawManager  An initialized MaterialLawManager object
        void init(Opm::DeckConstPtr deck,
                  std::shared_ptr<MaterialLawManager> materialLawManager)
        {
            init(Opm::phaseUsageFromDeck(deck), materialLawManager);
        }

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

        /// Update saturation state for the hysteresis tracking 
        /// \param[in]  n      Number of data points. 
        /// \param[in]  s      Array of nP saturation values.             
        void updateSatHyst(const int n,
                           const int* cells,
                           const double* s);

        /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
        /// \param[in]     cell  Cell index. 
        /// \param[in]     pcow  P_oil - P_water.
        /// \param[in/out] swat  Water saturation. / Possibly modified Water saturation.        
        void swatInitScaling(const int cell, 
                             const double pcow, 
                             double & swat);

        /// Returns a reference to the MaterialLawManager
        const MaterialLawManager& materialLawManager() const { return *materialLawManager_; }


    private:
        std::shared_ptr<MaterialLawManager> materialLawManager_;
        PhaseUsage phaseUsage_;
    };



} // namespace Opm

#endif // OPM_SATURATIONPROPSFROMDECK_HEADER_INCLUDED
