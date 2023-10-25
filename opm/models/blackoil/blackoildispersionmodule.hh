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
 * \brief Classes required for mechanical dispersion.
 */
#ifndef EWOMS_DISPERSION_MODULE_HH
#define EWOMS_DISPERSION_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <stdexcept>

namespace Opm {

/*!
 * \ingroup Dispersion
 * \class Opm::BlackOilDispersionModule
 * \brief Provides the auxiliary methods required for consideration of the
 * dispersion equation. 
 */
template <class TypeTag, bool enableDispersion>
class BlackOilDispersionModule;

template <class TypeTag, bool enableDispersion>
class BlackOilDispersionExtensiveQuantities;

/*!
 * \copydoc Opm::BlackOilDispersionModule
 */
template <class TypeTag>
class BlackOilDispersionModule<TypeTag, /*enableDispersion=*/false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    enum { numPhases = FluidSystem::numPhases };

public:
    using ExtensiveQuantities = BlackOilDispersionExtensiveQuantities<TypeTag,false>;
    /*!
     * \brief Adds the dispersive flux to the flux vector over a flux
     *        integration point.
     */
    template <class Context>
    static void addDispersiveFlux(RateVector&,
                                  const Context&,
                                  unsigned,
                                  unsigned)
    {}

    template<class FluidState, class Scalar>
    static void addDispersiveFlux(RateVector&,
                                  const FluidState&,
                                  const FluidState&,
                                  const Evaluation&,
                                  const Scalar&)
    {}
};

/*!
 * \copydoc Opm::BlackOilDispersionModule
 */
template <class TypeTag>
class BlackOilDispersionModule<TypeTag, /*enableDispersion=*/true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };

    using Toolbox = MathToolbox<Evaluation>;

public:
    using ExtensiveQuantities = BlackOilDispersionExtensiveQuantities<TypeTag,true>;
    /*!
     * \brief Adds the mass flux due to dispersion to the flux vector over the
     *        flux integration point.
     */
    template <class Context>
    static void addDispersiveFlux(RateVector& flux, const Context& context,
                                  unsigned spaceIdx, unsigned timeIdx)
    {
        // Only work if dispersion is enabled by DISPERC in the deck
        if (!context.simulator().vanguard().eclState().getSimulationConfig().rock_config().dispersion()) {
            return;
        }
        const auto& extQuants = context.extensiveQuantities(spaceIdx, timeIdx);
        const auto& fluidStateI = context.intensiveQuantities(extQuants.interiorIndex(), timeIdx).fluidState();
        const auto& fluidStateJ = context.intensiveQuantities(extQuants.exteriorIndex(), timeIdx).fluidState();
        const auto& dispersivity = extQuants.dispersivity();
        const auto& normVelocityAvg = extQuants.normVelocityAvg();
        addDispersiveFlux(flux, fluidStateI, fluidStateJ, dispersivity, normVelocityAvg);
    }

    /*!
     * \brief Adds the mass flux due to dispersion to the flux vector over the
     *        integration point. Following the notation in blackoilmodel.hh,
     *        the dispersive flux for component \f$\kappa\f$ in phase \f$\alpha\f$
     *        is given by: \f$-b_\alpha E||\mathrm{v}_\alpha||\mathbf{grad}X_\alpha^\kappa\f$,
     *        where \f$b_\alpha\f$ is the shrinkage/expansion factor [-], E is the isotropic
     *        dispersivity coefficient [L], \f$\mathrm{v}_\alpha\f$ is the filter velocity
     *        [L/T], and \f$X_\alpha^\kappa\f$ the component mass fraction [-]. Each component mass
     *        fraction can be computed using \f$R_s,\;R_v,\;R_{sw},\;R_{vw}\f$. For example,
     *        \f$X_w^G=\frac{R_{sw}}{R_{sw}+\rho_w/\rho_g}\f$, where \f$\rho_w\f$ and \f$\rho_g\f$
     *        are the reference densities.
     *        Following the implementation of the diffusive flux (blackoildiffusionmodule.hh) and considering
     *        the case for the water phase and gas component as an example, for cells i and j, the discrete version
     *        of the dispersive flux at the face's integration point is given by
     *        \f$-b_{w,ij}v_{w,ij}(\frac{1}{R_{sw,ij}+\rho_w/\rho_g})D_{ij}(R_{sw,i}-R_{sw,j})\f$
     *        where \f$b_{w,ij}\f$, \f$v_{w,ij}\f$, and \f$R_{sw,ij}\f$ are computed using the arithmetic mean, and
     *        the ratio \f$\frac{1}{R_{sw,ij}+\rho_w/\rho_g}\f$ is denoted as conversion factor. The dispersivity
     *        \f$D_{ij}\f$ is computed in ecltransmissibility_impl.hh, using the dispersion coefficients \f$E_i\f$
     *        and \f$E_j\f$.
     */
    template<class FluidState, class Scalar>
    static void addDispersiveFlux(RateVector& flux,
                                  const FluidState& fluidStateI,
                                  const FluidState& fluidStateJ,
                                  const Evaluation& dispersivity,
                                  const Scalar& normVelocityAvg)
    {
        unsigned pvtRegionIndex = fluidStateI.pvtRegionIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            // no dispersion in water for blackoil models unless water can contain dissolved gas
            if (!FluidSystem::enableDissolvedGasInWater() && FluidSystem::waterPhaseIdx == phaseIdx) {
                continue;
            }

            // no dispersion in gas for blackoil models unless gas can contain evaporated water or oil
            if ((!FluidSystem::enableVaporizedWater() && !FluidSystem::enableVaporizedOil()) && FluidSystem::gasPhaseIdx == phaseIdx) {
                continue;
            }

            // arithmetic mean of the phase's b factor
            Evaluation bAvg = fluidStateI.invB(phaseIdx);
            bAvg += Toolbox::value(fluidStateJ.invB(phaseIdx));
            bAvg /= 2;

            Evaluation convFactor = 1.0;
            Evaluation diffR = 0.0;
            if (FluidSystem::enableDissolvedGas() && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && phaseIdx == FluidSystem::oilPhaseIdx) {
                Evaluation rsAvg = (fluidStateI.Rs() + Toolbox::value(fluidStateJ.Rs())) / 2;
                convFactor = 1.0 / (toMassFractionGasOil(pvtRegionIndex) + rsAvg);
                diffR = fluidStateI.Rs() - Toolbox::value(fluidStateJ.Rs());
            }
            if (FluidSystem::enableVaporizedOil() && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && phaseIdx == FluidSystem::gasPhaseIdx) {
                Evaluation rvAvg = (fluidStateI.Rv() + Toolbox::value(fluidStateJ.Rv())) / 2;
                convFactor = toMassFractionGasOil(pvtRegionIndex) / (1.0 + rvAvg*toMassFractionGasOil(pvtRegionIndex));
                diffR = fluidStateI.Rv() - Toolbox::value(fluidStateJ.Rv());
            }
            if (FluidSystem::enableDissolvedGasInWater() && phaseIdx == FluidSystem::waterPhaseIdx) {
                Evaluation rsAvg = (fluidStateI.Rsw() + Toolbox::value(fluidStateJ.Rsw())) / 2;
                convFactor = 1.0 / (toMassFractionGasWater(pvtRegionIndex) + rsAvg);
                diffR = fluidStateI.Rsw() - Toolbox::value(fluidStateJ.Rsw());
            }
            if (FluidSystem::enableVaporizedWater() && phaseIdx == FluidSystem::gasPhaseIdx) {
                Evaluation rvAvg = (fluidStateI.Rvw() + Toolbox::value(fluidStateJ.Rvw())) / 2;
                convFactor = toMassFractionGasWater(pvtRegionIndex)/ (1.0 + rvAvg*toMassFractionGasWater(pvtRegionIndex));
                diffR = fluidStateI.Rvw() - Toolbox::value(fluidStateJ.Rvw());
            }

            // mass flux of solvent component
            unsigned solventCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            unsigned activeSolventCompIdx = Indices::canonicalToActiveComponentIndex(solventCompIdx);
            flux[conti0EqIdx + activeSolventCompIdx] +=
                    - bAvg
                    * normVelocityAvg[phaseIdx]
                    * convFactor
                    * dispersivity
                    * diffR;

            // mass flux of solute component
            unsigned soluteCompIdx = FluidSystem::soluteComponentIndex(phaseIdx);
            unsigned activeSoluteCompIdx = Indices::canonicalToActiveComponentIndex(soluteCompIdx);
            flux[conti0EqIdx + activeSoluteCompIdx] +=
                    bAvg
                    * normVelocityAvg[phaseIdx]
                    * convFactor
                    * dispersivity
                    * diffR;
        }
    }

private:
    static Scalar toMassFractionGasOil (unsigned regionIdx) {
        Scalar rhoO = FluidSystem::referenceDensity(FluidSystem::oilPhaseIdx, regionIdx);
        Scalar rhoG = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);
        return rhoO / rhoG;
    }
    static Scalar toMassFractionGasWater (unsigned regionIdx) {
        Scalar rhoW = FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, regionIdx);
        Scalar rhoG = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);
        return rhoW / rhoG;
    }
};

/*!
 * \ingroup Dispersion
 * \class Opm::BlackOilDispersionIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the
 *        calculation of dispersive fluxes.
 */
template <class TypeTag, bool enableDispersion>
class BlackOilDispersionIntensiveQuantities;

/*!
 * \copydoc Opm::DispersionIntensiveQuantities
 */
template <class TypeTag>
class BlackOilDispersionIntensiveQuantities<TypeTag, /*enableDispersion=*/false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    /*!
     * \brief Returns the max. norm of the filter velocity of the cell.
     */
    Scalar normVelocityCell(unsigned, unsigned) const
    {
        throw std::logic_error("Method normVelocityCell() "
                               "does not make sense if dispersion is disabled");
    }

protected:
    /*!
     * \brief Update the quantities required to calculate dispersive
     *        fluxes.
     */
    template<class ElementContext>
    void update_(ElementContext&,
                 unsigned,
                 unsigned)
    { }
};

/*!
 * \copydoc Opm::DispersionIntensiveQuantities
 */
template <class TypeTag>
class BlackOilDispersionIntensiveQuantities<TypeTag, /*enableDispersion=*/true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };

public:    
    /*!
     * \brief Returns the max. norm of the filter velocity of the cell.
     */
    Scalar normVelocityCell(unsigned phaseIdx) const
    {

        return normVelocityCell_[phaseIdx];
    }

protected:
    /*!
     * \brief Update the quantities required to calculate dispersive
     *        mass fluxes. This considers the linear disperison model
     *        described in the SPE CSP11 benchmark document (eq. 2.3)
     *        https://github.com/Simulation-Benchmarks/11thSPE-CSP/
     *        blob/main/description/spe_csp11_description.pdf
     *        The maximum norm is used to compute the cell
     *        filter velocity value of the corresponding phase.
     */
    template<class ElementContext>
    void update_(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        // Only work if dispersion is enabled by DISPERC in the deck
        if (!elemCtx.simulator().vanguard().eclState().getSimulationConfig().rock_config().dispersion()) {
            return;
        }
        const auto& problem = elemCtx.simulator().problem();
        if (problem.model().linearizer().getVelocityInfo().empty()) {
            return;
        }
        const std::array<int, 3> phaseIdxs = { gasPhaseIdx, oilPhaseIdx, waterPhaseIdx };
        const std::array<int, 3> compIdxs = { gasCompIdx, oilCompIdx, waterCompIdx };
        const auto& velocityInf = problem.model().linearizer().getVelocityInfo();
        unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        auto velocityInfos = velocityInf[globalDofIdx];
        for (unsigned i = 0; i < phaseIdxs.size(); ++i) {
            normVelocityCell_[i] = 0;
        }
        for (auto& velocityInfo : velocityInfos) {
            for (unsigned i = 0; i < phaseIdxs.size(); ++i) {
                if (FluidSystem::phaseIsActive(phaseIdxs[i])) {
                    normVelocityCell_[phaseIdxs[i]] = max( normVelocityCell_[phaseIdxs[i]], 
                        std::abs( velocityInfo.velocity[conti0EqIdx
                            + Indices::canonicalToActiveComponentIndex(compIdxs[i])] ));
                }
            }
        }
    }

private:
    Scalar normVelocityCell_[numPhases];
};

/*!
 * \ingroup Dispersion
 * \class Opm::BlackOilDispersionExtensiveQuantities
 *
 * \brief Provides the quantities required to calculate dispersive mass fluxes.
 */
template <class TypeTag, bool enableDispersion>
class BlackOilDispersionExtensiveQuantities;

/*!
 * \copydoc Opm::DispersionExtensiveQuantities
 */
template <class TypeTag>
class BlackOilDispersionExtensiveQuantities<TypeTag, /*enableDispersion=*/false>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    enum { numPhases = FluidSystem::numPhases };

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        the dispersive fluxes.
     */
    void update_(const ElementContext&,
                 unsigned,
                 unsigned)
    {}

    template <class Context, class FluidState>
    void updateBoundary_(const Context&,
                         unsigned,
                         unsigned,
                         const FluidState&)
    {}

public:
    using ScalarArray = Scalar[numPhases];

    static void update(ScalarArray&,
                       const IntensiveQuantities&,
                       const IntensiveQuantities&)
    {}

    /*!
     * \brief The dispersivity the face.
     *
     */
    const Scalar& dispersivity() const
    {
        throw std::logic_error("The method dispersivity() does not "
                               "make sense if dispersion is disabled.");
    }

    /*!
     * \brief The effective filter velocity coefficient in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Scalar& normVelocityAvg(unsigned) const
    {
        throw std::logic_error("The method normVelocityAvg() "
                               "does not make sense if dispersion is disabled.");
    }

};

/*!
 * \copydoc Opm::BlackOilDispersionExtensiveQuantities
 */
template <class TypeTag>
class BlackOilDispersionExtensiveQuantities<TypeTag, /*enableDispersion=*/true>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Toolbox = MathToolbox<Evaluation>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;
public:
    using ScalarArray = Scalar[numPhases];
    static void update(ScalarArray& normVelocityAvg,
                       const IntensiveQuantities& intQuantsInside,
                       const IntensiveQuantities& intQuantsOutside)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            // no dispersion in water for blackoil models unless water can contain dissolved gas
            if (!FluidSystem::enableDissolvedGasInWater() && FluidSystem::waterPhaseIdx == phaseIdx) {
                continue;
            }
            // no dispersion in gas for blackoil models unless gas can contain evaporated water or oil
            if ((!FluidSystem::enableVaporizedWater() && !FluidSystem::enableVaporizedOil()) && FluidSystem::gasPhaseIdx == phaseIdx) {
                continue;
            }
            // use the arithmetic average for the effective
            // velocity coefficients at the face's integration point.
            normVelocityAvg[phaseIdx] = 0.5 *
                ( intQuantsInside.normVelocityCell(phaseIdx) +
                    intQuantsOutside.normVelocityCell(phaseIdx) );
            Valgrind::CheckDefined(normVelocityAvg[phaseIdx]);
        }
    }
protected:
    template <class Context, class FluidState>
    void updateBoundary_(const Context&,
                         unsigned,
                         unsigned,
                         const FluidState&)
    {
        throw std::runtime_error("Not implemented: Dispersion across boundary not implemented for blackoil");
    }

public:
    /*!
     * \brief The dispersivity of the face.
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Scalar& dispersivity() const
    { return dispersivity_; }

    /*!
     * \brief The effective velocity coefficient in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Scalar& normVelocityAvg(unsigned phaseIdx) const
    { return normVelocityAvg_[phaseIdx]; }

    const auto& normVelocityAvg() const{
        return normVelocityAvg_;
    }

private:
    Scalar dispersivity_;
    ScalarArray normVelocityAvg_;
};

} // namespace Opm

#endif
