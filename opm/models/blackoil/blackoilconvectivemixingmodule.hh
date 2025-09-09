// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Classes required for dynamic convective mixing.
 */
#ifndef EWOMS_CONVECTIVEMIXING_MODULE_HH
#define EWOMS_CONVECTIVEMIXING_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#if HAVE_ECL_INPUT
#include <cstddef>
#endif

#include <vector>

namespace Opm {

/*!
 * \copydoc Opm::BlackOilConvectiveMixingModule
 * \brief Provides the convective term in the transport flux for the brine
 * when convective mixing (enhanced dissolution of CO2 in brine) occurs.
 * Controlled by the regimes for a controlvolume:
 * i) initial phase (CO2 dissolves in brine due to diffusion)
 * ii) linear phase (Convective fingers of CO2-rich brine propagate downwards)
 * iii) steady-state-phase (fingers have passed through the bottom of a control
 * -volume but the larger scale convective process is still active)
 * iv) decline phase (Convection ceases at the large-scale when the CO2
 * has been completely dissolved)
 */

template <class TypeTag, bool enableConvectiveMixing>
class BlackOilConvectiveMixingModule;

/*!
 * \copydoc Opm::BlackOilConvectiveMixingModule
 */

template <class TypeTag>
class BlackOilConvectiveMixingModule<TypeTag, /*enableConvectiveMixing=*/false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { dimWorld = GridView::dimensionworld };

public:
    struct ConvectiveMixingModuleParam
    {};

    #if HAVE_ECL_INPUT
    static void beginEpisode(const EclipseState&,
                             const Schedule&,
                             const int,
                             ConvectiveMixingModuleParam&)
    {}
    #endif

    template <class Context>
    static bool active(const Context&)
    { return false; }

    static void modifyAvgDensity(Evaluation&,
                                 const IntensiveQuantities&,
                                 const IntensiveQuantities&,
                                 const unsigned int,
                                 const ConvectiveMixingModuleParam&)
    {}

    template <class Context>
    static void addConvectiveMixingFlux(RateVector&,
                                        const Context&,
                                        unsigned,
                                        unsigned)
    {}

    /*!
     * \brief Adds the convective mixing mass flux flux to the flux vector over a flux
     *        integration point.
     */
    static void addConvectiveMixingFlux(RateVector&,
                                        const IntensiveQuantities&,
                                        const IntensiveQuantities&,
                                        const unsigned,
                                        const unsigned,
                                        const Scalar,
                                        const Scalar,
                                        const Scalar,
                                        const ConvectiveMixingModuleParam&)
    {}
};

template <class TypeTag>
class BlackOilConvectiveMixingModule<TypeTag, /*enableConvectiveMixing=*/true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Toolbox = MathToolbox<Evaluation>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { dimWorld = GridView::dimensionworld };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    static constexpr unsigned contiEnergyEqIdx = Indices::contiEnergyEqIdx;

public:
    struct ConvectiveMixingModuleParam
    {
        std::vector<bool> active_;
        std::vector<Scalar> Xhi_;
        std::vector<Scalar> Psi_;
    };

    #if HAVE_ECL_INPUT
    static void beginEpisode(const EclipseState& eclState,
                             const Schedule& schedule,
                             const int episodeIdx,
                             ConvectiveMixingModuleParam& info)
    {
        // check that Xhi and Psi didn't change
        std::size_t numRegions = eclState.runspec().tabdims().getNumPVTTables();
        const auto& control = schedule[episodeIdx].oilvap();
        if (info.active_.empty()) {
            info.active_.resize(numRegions);
            info.Xhi_.resize(numRegions);
            info.Psi_.resize(numRegions);
        }
        for (std::size_t i = 0; i < numRegions; ++i ) {
            info.active_[i] = control.drsdtConvective(i);
            if (control.drsdtConvective(i)) {
                info.Xhi_[i] = control.getMaxDRSDT(i);
                info.Psi_[i] = control.getPsi(i);
            }
        }
    }
    #endif

    static void modifyAvgDensity(Evaluation& rhoAvg,
                                 const IntensiveQuantities& intQuantsIn,
                                 const IntensiveQuantities& intQuantsEx,
                                 const unsigned phaseIdx,
                                 const ConvectiveMixingModuleParam& info) {
        
        if (info.active_.empty()) {
            return;
        }
        if (!info.active_[intQuantsIn.pvtRegionIndex()] || !info.active_[intQuantsEx.pvtRegionIndex()]) {
            return;
        }

        if (phaseIdx == FluidSystem::gasPhaseIdx) {
            return;
        }

        const auto& liquidPhaseIdx =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
               ? FluidSystem::waterPhaseIdx
               : FluidSystem::oilPhaseIdx;

        // Compute avg density based on pure water
        const auto& t_in = intQuantsIn.fluidState().temperature(liquidPhaseIdx);
        const auto& p_in = intQuantsIn.fluidState().pressure(liquidPhaseIdx);
        const auto& salt_in = intQuantsIn.fluidState().saltConcentration();

        const auto& bLiquidIn =
            FluidSystem::phaseIsActive(waterPhaseIdx)
                ? FluidSystem::waterPvt().inverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(),
                                                                         t_in, p_in, Evaluation(0.0), salt_in)
                : FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(),
                                                                 t_in, p_in, Evaluation(0.0));

        const auto& refDensityLiquidIn =
            FluidSystem::phaseIsActive(waterPhaseIdx)
                ? FluidSystem::waterPvt().waterReferenceDensity(intQuantsIn.pvtRegionIndex())
                : FluidSystem::oilPvt().oilReferenceDensity(intQuantsIn.pvtRegionIndex());

        const auto& rho_in =  bLiquidIn * refDensityLiquidIn;

        const auto t_ex = Toolbox::value(intQuantsEx.fluidState().temperature(liquidPhaseIdx));
        const auto p_ex = Toolbox::value(intQuantsEx.fluidState().pressure(liquidPhaseIdx));
        const auto salt_ex = Toolbox::value(intQuantsEx.fluidState().saltConcentration());

        const auto bLiquidEx =
            FluidSystem::phaseIsActive(waterPhaseIdx)
                ? FluidSystem::waterPvt().inverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(),
                                                                       t_ex, p_ex, Scalar{0.0}, salt_ex)
                : FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(),
                                                                     t_ex, p_ex, Scalar{0.0});

        const auto& refDensityLiquidEx =
            FluidSystem::phaseIsActive(waterPhaseIdx)
                ? FluidSystem::waterPvt().waterReferenceDensity(intQuantsEx.pvtRegionIndex())
                : FluidSystem::oilPvt().oilReferenceDensity(intQuantsEx.pvtRegionIndex());

        const auto rho_ex =  bLiquidEx * refDensityLiquidEx;

        rhoAvg = (rho_in + rho_ex) / 2;
    }

    template <class Context>
    static void addConvectiveMixingFlux(RateVector& flux,
                                        const Context& elemCtx,
                                        unsigned scvfIdx,
                                        unsigned timeIdx)
    {
        // need for darcy flux calculation
        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        unsigned interiorDofIdx = scvf.interiorIndex();
        unsigned exteriorDofIdx = scvf.exteriorIndex();
        assert(interiorDofIdx != exteriorDofIdx);

        const auto& globalIndexIn = stencil.globalSpaceIndex(interiorDofIdx);
        const auto& globalIndexEx = stencil.globalSpaceIndex(exteriorDofIdx);
        Scalar trans = problem.transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar faceArea = scvf.area();
        const Scalar g = problem.gravity()[dimWorld - 1];
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);
        const Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        const Scalar zEx = problem.dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);
        const Scalar distZ = zIn - zEx;
        addConvectiveMixingFlux(flux,
                                intQuantsIn,
                                intQuantsEx,
                                globalIndexIn,
                                globalIndexEx,
                                distZ * g,
                                trans,
                                faceArea,
                                problem.moduleParams().convectiveMixingModuleParam);
    }

    /*!
     * \brief Adds the convective mixing mass flux flux to the flux vector over a flux
     *        integration point.
      */
    static void addConvectiveMixingFlux(RateVector& flux,
                                        const IntensiveQuantities& intQuantsIn,
                                        const IntensiveQuantities& intQuantsEx,
                                        const unsigned globalIndexIn,
                                        const unsigned globalIndexEx,
                                        const Scalar distZg,
                                        const Scalar trans,
                                        const Scalar faceArea,
                                        const ConvectiveMixingModuleParam& info)
    {
        if (info.active_.empty()) {
            return;
        }

        if (!info.active_[intQuantsIn.pvtRegionIndex()] || !info.active_[intQuantsEx.pvtRegionIndex()]) {
            return;
        }

        const auto& liquidPhaseIdx =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPhaseIdx
                : FluidSystem::oilPhaseIdx;

        // interiour
        const auto& t_in = intQuantsIn.fluidState().temperature(liquidPhaseIdx);
        const auto& p_in = intQuantsIn.fluidState().pressure(liquidPhaseIdx);
        const auto& rssat_in = intQuantsIn.saturatedDissolutionFactor();
        const auto& salt_in = intQuantsIn.fluidState().saltSaturation();

        const auto bLiquidSatIn =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPvt().inverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(),
                                                                       t_in, p_in, rssat_in, salt_in)
                : FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(),
                                                                     t_in, p_in, rssat_in);

        const auto& densityLiquidIn =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPvt().waterReferenceDensity(intQuantsIn.pvtRegionIndex())
                : FluidSystem::oilPvt().oilReferenceDensity(intQuantsIn.pvtRegionIndex());

        const auto rho_in = Opm::getValue(intQuantsIn.fluidState().invB(liquidPhaseIdx)) * densityLiquidIn;
        const auto rho_sat_in = bLiquidSatIn *
                                (densityLiquidIn +
                                 rssat_in * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx,
                                                                          intQuantsIn.pvtRegionIndex()));

        // exteriour
        const auto t_ex = Opm::getValue(intQuantsEx.fluidState().temperature(liquidPhaseIdx));
        const auto p_ex = Opm::getValue(intQuantsEx.fluidState().pressure(liquidPhaseIdx));
        const auto rssat_ex = Opm::getValue(intQuantsEx.saturatedDissolutionFactor());
        const auto salt_ex = Opm::getValue(intQuantsEx.fluidState().saltSaturation());
        const auto bLiquidSatEx =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPvt().inverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(),
                                                                       t_ex, p_ex, rssat_ex, salt_ex)
                : FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(),
                                                                     t_ex, p_ex, rssat_ex);

        const auto& densityLiquidEx =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPvt().waterReferenceDensity(intQuantsEx.pvtRegionIndex())
                : FluidSystem::oilPvt().oilReferenceDensity(intQuantsEx.pvtRegionIndex());

        const auto rho_ex = Opm::getValue(intQuantsEx.fluidState().invB(liquidPhaseIdx)) * densityLiquidEx;
        const auto rho_sat_ex = bLiquidSatEx *
                                (densityLiquidEx +
                                 rssat_ex * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx,
                                                                          intQuantsEx.pvtRegionIndex()));
        // rho difference approximation
        const auto delta_rho = (rho_sat_ex + rho_sat_in - rho_in - rho_ex) / 2;
        const auto pressure_difference_convective_mixing =  delta_rho * distZg;

        // if change in pressure
        if (Opm::abs(pressure_difference_convective_mixing) > 1e-12) {
            // find new upstream direction
            short interiorDofIdx = 0;
            constexpr short exteriorDofIdx = 1;
            short upIdx = 0;

            if (pressure_difference_convective_mixing > 0) {
                upIdx = exteriorDofIdx;
            }

            const auto& up = (upIdx == interiorDofIdx) ? intQuantsIn : intQuantsEx;
            const auto& rssat_up = (upIdx == interiorDofIdx) ? rssat_in : rssat_ex;
            unsigned globalUpIndex = (upIdx == interiorDofIdx) ? globalIndexIn : globalIndexEx;
            const auto& Rsup =
                FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                    ? up.fluidState().Rsw()
                    : up.fluidState().Rs();

            const Evaluation& transMult = up.rockCompTransMultiplier();
            const auto& invB = up.fluidState().invB(liquidPhaseIdx);
            const auto& visc = up.fluidState().viscosity(liquidPhaseIdx);

            // We restrict the convective mixing mass flux to rssat * Psi.
            const Evaluation RsupRestricted = Opm::min(Rsup, rssat_up*info.Psi_[up.pvtRegionIndex()]);

            const auto convectiveFlux = -trans * transMult * info.Xhi_[up.pvtRegionIndex()] * invB *
                                        pressure_difference_convective_mixing * RsupRestricted / (visc * faceArea);
            unsigned activeGasCompIdx = FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx);
            if (globalUpIndex == globalIndexIn) {
                flux[conti0EqIdx + activeGasCompIdx] += convectiveFlux;
            }
            else {
                flux[conti0EqIdx + activeGasCompIdx] += Opm::getValue(convectiveFlux);
            }

            if constexpr (enableEnergy) {
                const auto& h = up.fluidState().enthalpy(liquidPhaseIdx) *
                                FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, up.pvtRegionIndex());
                if (globalUpIndex == globalIndexIn) {
                    flux[contiEnergyEqIdx] += convectiveFlux * h;
                }
                else {
                    flux[contiEnergyEqIdx] += Opm::getValue(h) * Opm::getValue(convectiveFlux);
                }
            }
        }
    }
};

template <class TypeTag, bool enableConvectiveMixingV>
class BlackOilConvectiveMixingIntensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilConvectiveMixingIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        convective mixing (DRSDTCON) model.
 */
template <class TypeTag>
class BlackOilConvectiveMixingIntensiveQuantities<TypeTag, /*enableConvectiveMixingV=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    /*!
     * \brief Compute the intensive quantities needed to handle convective dissolution
     *
     */
    void updateSaturatedDissolutionFactor_()
    {
        const auto liquidPhaseIdx =
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? FluidSystem::waterPhaseIdx
                : FluidSystem::oilPhaseIdx;

        const Evaluation SoMax = 0.0;
        saturatedDissolutionFactor_ = FluidSystem::saturatedDissolutionFactor(asImp_().fluidState(),
                                                                              liquidPhaseIdx,
                                                                              asImp_().pvtRegionIndex(),
                                                                              SoMax);
    }

    const Evaluation& saturatedDissolutionFactor() const
    { return saturatedDissolutionFactor_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation saturatedDissolutionFactor_;
};

template <class TypeTag>
class BlackOilConvectiveMixingIntensiveQuantities<TypeTag, false>
{
};

}

#endif
