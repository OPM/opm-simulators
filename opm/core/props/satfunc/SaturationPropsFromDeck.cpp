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
#include "config.h"

#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <opm/core/utility/UniformTableLinear.hpp>
#include <opm/core/utility/NonuniformTableLinear.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/simulator/ExplicitArraysFluidState.hpp>
#include <opm/core/simulator/ExplicitArraysSatDerivativesFluidState.hpp>

#include <iostream>
#include <map>

namespace Opm
{

    typedef SaturationPropsFromDeck::MaterialLawManager::MaterialLaw MaterialLaw;

    // ----------- Methods of SaturationPropsFromDeck ---------


    /// Default constructor.
    SaturationPropsFromDeck::SaturationPropsFromDeck()
    {
    }

    /// Initialize from deck.
    void SaturationPropsFromDeck::init(const PhaseUsage &phaseUsage,
                                       std::shared_ptr<MaterialLawManager> materialLawManager)
    {
        phaseUsage_ = phaseUsage;
        materialLawManager_ = materialLawManager;
    }

    /// \return   P, the number of phases.
    int SaturationPropsFromDeck::numPhases() const
    {
        return phaseUsage_.num_phases;
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
    void SaturationPropsFromDeck::relperm(const int n,
                                          const double* s,
                                          const int* cells,
                                          double* kr,
                                          double* dkrds) const
    {
        assert(cells != 0);

        const int np = numPhases();
        if (dkrds) {
            ExplicitArraysSatDerivativesFluidState fluidState(phaseUsage_);
            fluidState.setSaturationArray(s);

            typedef ExplicitArraysSatDerivativesFluidState::Evaluation Evaluation;
            Evaluation relativePerms[BlackoilPhases::MaxNumPhases];
            for (int i = 0; i < n; ++i) {
                fluidState.setIndex(i);
                const auto& params = materialLawManager_->materialLawParams(cells[i]);
                MaterialLaw::relativePermeabilities(relativePerms, params, fluidState);

                // copy the values calculated using opm-material to the target arrays
                for (int krPhaseIdx = 0; krPhaseIdx < np; ++krPhaseIdx) {
                    kr[np*i + krPhaseIdx] = relativePerms[krPhaseIdx].value;

                    for (int satPhaseIdx = 0; satPhaseIdx < np; ++satPhaseIdx)
                        dkrds[np*np*i + satPhaseIdx*np + krPhaseIdx] = relativePerms[krPhaseIdx].derivatives[satPhaseIdx];
                }
            }
        } else {
            ExplicitArraysFluidState fluidState(phaseUsage_);
            fluidState.setSaturationArray(s);

            double relativePerms[BlackoilPhases::MaxNumPhases];
            for (int i = 0; i < n; ++i) {
                fluidState.setIndex(i);
                const auto& params = materialLawManager_->materialLawParams(cells[i]);
                MaterialLaw::relativePermeabilities(relativePerms, params, fluidState);

                // copy the values calculated using opm-material to the target arrays
                for (int krPhaseIdx = 0; krPhaseIdx < np; ++krPhaseIdx) {
                    kr[np*i + krPhaseIdx] = relativePerms[krPhaseIdx];
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
    void SaturationPropsFromDeck::capPress(const int n,
                                           const double* s,
                                           const int* cells,
                                           double* pc,
                                           double* dpcds) const
    {
        assert(cells != 0);        
        assert(phaseUsage_.phase_used[BlackoilPhases::Liquid]);

        const int np = numPhases();

        if (dpcds) {
            ExplicitArraysSatDerivativesFluidState fluidState(phaseUsage_);
            typedef ExplicitArraysSatDerivativesFluidState::Evaluation Evaluation;
            fluidState.setSaturationArray(s);

            Evaluation capillaryPressures[BlackoilPhases::MaxNumPhases];
            for (int i = 0; i < n; ++i) {
                fluidState.setIndex(i);
                const auto& params = materialLawManager_->materialLawParams(cells[i]);
                MaterialLaw::capillaryPressures(capillaryPressures, params, fluidState);

                // copy the values calculated using opm-material to the target arrays
                for (int canonicalPhaseIdx = 0; canonicalPhaseIdx < BlackoilPhases::MaxNumPhases; ++canonicalPhaseIdx) {
                    // skip unused phases
                    if ( ! phaseUsage_.phase_used[canonicalPhaseIdx]) {
                        continue;
                    }
                    const int pcPhaseIdx = phaseUsage_.phase_pos[canonicalPhaseIdx];

                    const double sign = (canonicalPhaseIdx == BlackoilPhases::Aqua)? -1.0 : 1.0;
                    // in opm-material the wetting phase is the reference phase
                    // for two-phase problems i.e water for oil-water system,
                    // but for flow it is always oil. Add oil (liquid) capillary pressure value
                    // to shift the reference phase to oil
                    pc[np*i + pcPhaseIdx] = capillaryPressures[BlackoilPhases::Liquid].value + sign * capillaryPressures[canonicalPhaseIdx].value;
                    for (int canonicalSatPhaseIdx = 0; canonicalSatPhaseIdx < BlackoilPhases::MaxNumPhases; ++canonicalSatPhaseIdx) {
                        if ( ! phaseUsage_.phase_used[canonicalSatPhaseIdx])
                            continue;

                        const int satPhaseIdx = phaseUsage_.phase_pos[canonicalSatPhaseIdx];
                        dpcds[np*np*i + satPhaseIdx*np + pcPhaseIdx] = capillaryPressures[BlackoilPhases::Liquid].derivatives[canonicalSatPhaseIdx] + sign * capillaryPressures[canonicalPhaseIdx].derivatives[canonicalSatPhaseIdx];
                    }
                }
            }
        } else {
            ExplicitArraysFluidState fluidState(phaseUsage_);
            fluidState.setSaturationArray(s);

            double capillaryPressures[BlackoilPhases::MaxNumPhases];
            for (int i = 0; i < n; ++i) {         
                fluidState.setIndex(i);
                const auto& params = materialLawManager_->materialLawParams(cells[i]);
                MaterialLaw::capillaryPressures(capillaryPressures, params, fluidState);

                // copy the values calculated using opm-material to the target arrays
                for (int canonicalPhaseIdx = 0; canonicalPhaseIdx < BlackoilPhases::MaxNumPhases; ++canonicalPhaseIdx) {
                    // skip unused phases
                    if ( ! phaseUsage_.phase_used[canonicalPhaseIdx])
                        continue;

                    const int pcPhaseIdx = phaseUsage_.phase_pos[canonicalPhaseIdx];
                    double sign = (canonicalPhaseIdx == BlackoilPhases::Aqua)? -1.0 : 1.0;
                    // in opm-material the wetting phase is the reference phase
                    // for two-phase problems i.e water for oil-water system,
                    // but for flow it is always oil. Add oil (liquid) capillary pressure value
                    // to shift the reference phase to oil
                    pc[np*i + pcPhaseIdx] = capillaryPressures[BlackoilPhases::Liquid] + sign * capillaryPressures[canonicalPhaseIdx];
                }
            }
        }
    }


    /// Obtain the range of allowable saturation values.
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void SaturationPropsFromDeck::satRange(const int n,
                                           const int* cells,
                                           double* smin,
                                           double* smax) const
    {
        int wpos = phaseUsage_.phase_pos[BlackoilPhases::Aqua];
        int gpos = phaseUsage_.phase_pos[BlackoilPhases::Vapour];
        int opos = phaseUsage_.phase_pos[BlackoilPhases::Liquid];

        const int np = numPhases();
        for (int i = 0; i < n; ++i) {
            const auto& scaledDrainageInfo =
                materialLawManager_->oilWaterScaledEpsInfoDrainage(cells[i]);

            if (phaseUsage_.phase_used[BlackoilPhases::Aqua]) {
                smin[np*i + wpos] = scaledDrainageInfo.Swl;
                smax[np*i + wpos] = scaledDrainageInfo.Swu;
            }

            if (phaseUsage_.phase_used[BlackoilPhases::Vapour]) {
                smin[np*i + gpos] = scaledDrainageInfo.Sgl;
                smax[np*i + gpos] = scaledDrainageInfo.Sgu;
            }

            if (phaseUsage_.phase_used[BlackoilPhases::Liquid]) {
                smin[np*i + opos] = 1.0;
                smax[np*i + opos] = 1.0;
                if (phaseUsage_.phase_used[BlackoilPhases::Aqua]) {
                    smin[np*i + opos] -= smax[np*i + wpos];
                    smax[np*i + opos] -= smin[np*i + wpos];
                }
                if (phaseUsage_.phase_used[BlackoilPhases::Vapour]) {
                    smin[np*i + opos] -= smax[np*i + gpos];
                    smax[np*i + opos] -= smin[np*i + gpos];
                }
                smin[np*i + opos] = std::max(0.0, smin[np*i + opos]);
            }
        }
    }

    /// Update saturation state for the hysteresis tracking 
    /// \param[in]  n      Number of data points. 
    /// \param[in]  s      Array of nP saturation values.
    void SaturationPropsFromDeck::updateSatHyst(const int n,
                                                            const int* cells,
                                                            const double* s)
    {        
        assert(cells != 0);

        if (materialLawManager_->enableHysteresis()) {
            ExplicitArraysFluidState fluidState(phaseUsage_);
            fluidState.setSaturationArray(s);

            for (int i = 0; i < n; ++i) {
                fluidState.setIndex(i);
                materialLawManager_->updateHysteresis(fluidState, cells[i]);
            }
        }
    }


    /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
    /// \param[in]     cell  Cell index.
    /// \param[in]     pcow  P_oil - P_water.
    /// \param[in/out] swat  Water saturation. / Possibly modified Water saturation.
    void SaturationPropsFromDeck::swatInitScaling(const int cell,
                                                              const double pcow,
                                                              double& swat)
    {
        swat = materialLawManager_->applySwatinit(cell, pcow, swat);
    }
} // namespace Opm
