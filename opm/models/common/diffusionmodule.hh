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
 * \brief Classes required for molecular diffusion.
 */
#ifndef EWOMS_DIFFUSION_MODULE_HH
#define EWOMS_DIFFUSION_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(Indices);

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup Diffusion
 * \class Opm::DiffusionModule
 * \brief Provides the auxiliary methods required for consideration of the
 * diffusion equation.
 */
template <class TypeTag, bool enableDiffusion>
class DiffusionModule;

/*!
 * \copydoc Opm::DiffusionModule
 */
template <class TypeTag>
class DiffusionModule<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

public:
    /*!
     * \brief Register all run-time parameters for the diffusion module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Adds the diffusive mass flux flux to the flux vector over a flux
     *        integration point.
      */
    template <class Context>
    static void addDiffusiveFlux(RateVector& flux OPM_UNUSED,
                                 const Context& context OPM_UNUSED,
                                 unsigned spaceIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED)
    {}
};

/*!
 * \copydoc Opm::DiffusionModule
 */
template <class TypeTag>
class DiffusionModule<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \brief Register all run-time parameters for the diffusion module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Adds the mass flux due to molecular diffusion to the flux vector over the
     *        flux integration point.
     */
    template <class Context>
    static void addDiffusiveFlux(RateVector& flux, const Context& context,
                                 unsigned spaceIdx, unsigned timeIdx)
    {
        const auto& extQuants = context.extensiveQuantities(spaceIdx, timeIdx);

        const auto& fluidStateI = context.intensiveQuantities(extQuants.interiorIndex(), timeIdx).fluidState();
        const auto& fluidStateJ = context.intensiveQuantities(extQuants.exteriorIndex(), timeIdx).fluidState();

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // arithmetic mean of the phase's molar density
            Evaluation rhoMolar = fluidStateI.molarDensity(phaseIdx);
            rhoMolar += Toolbox::value(fluidStateJ.molarDensity(phaseIdx));
            rhoMolar /= 2;

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                // mass flux due to molecular diffusion
                flux[conti0EqIdx + compIdx] +=
                    -rhoMolar
                    * extQuants.moleFractionGradientNormal(phaseIdx, compIdx)
                    * extQuants.effectiveDiffusionCoefficient(phaseIdx, compIdx);
        }
    }
};

/*!
 * \ingroup Diffusion
 * \class Opm::DiffusionIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the
 *        calculation of molecular diffusive fluxes.
 */
template <class TypeTag, bool enableDiffusion>
class DiffusionIntensiveQuantities;

/*!
 * \copydoc Opm::DiffusionIntensiveQuantities
 */
template <class TypeTag>
class DiffusionIntensiveQuantities<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief Returns the tortuousity of the sub-domain of a fluid
     *        phase in the porous medium.
     */
    Scalar tortuosity(unsigned phaseIdx OPM_UNUSED) const
    {
        throw std::logic_error("Method tortuosity() does not make sense "
                               "if diffusion is disabled");
    }

    /*!
     * \brief Returns the molecular diffusion coefficient for a
     *        component in a phase.
     */
    Scalar diffusionCoefficient(unsigned phaseIdx OPM_UNUSED, unsigned compIdx OPM_UNUSED) const
    {
        throw std::logic_error("Method diffusionCoefficient() does not "
                               "make sense if diffusion is disabled");
    }

    /*!
     * \brief Returns the effective molecular diffusion coefficient of
     *        the porous medium for a component in a phase.
     */
    Scalar effectiveDiffusionCoefficient(unsigned phaseIdx OPM_UNUSED, unsigned compIdx OPM_UNUSED) const
    {
        throw std::logic_error("Method effectiveDiffusionCoefficient() "
                               "does not make sense if diffusion is disabled");
    }

protected:
    /*!
     * \brief Update the quantities required to calculate diffusive
     *        mass fluxes.
     */
    template <class FluidState>
    void update_(FluidState& fs OPM_UNUSED,
                 typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache OPM_UNUSED,
                 const ElementContext& elemCtx OPM_UNUSED,
                 unsigned dofIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED)
    { }
};

/*!
 * \copydoc Opm::DiffusionIntensiveQuantities
 */
template <class TypeTag>
class DiffusionIntensiveQuantities<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

public:
    /*!
     * \brief Returns the molecular diffusion coefficient for a
     *        component in a phase.
     */
    Evaluation diffusionCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return diffusionCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the tortuousity of the sub-domain of a fluid
     *        phase in the porous medium.
     */
    Evaluation tortuosity(unsigned phaseIdx) const
    { return tortuosity_[phaseIdx]; }

    /*!
     * \brief Returns the effective molecular diffusion coefficient of
     *        the porous medium for a component in a phase.
     */
    Evaluation effectiveDiffusionCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return tortuosity_[phaseIdx] * diffusionCoefficient_[phaseIdx][compIdx]; }

protected:
    /*!
     * \brief Update the quantities required to calculate diffusive
     *        mass fluxes.
     */
    template <class FluidState>
    void update_(FluidState& fluidState,
                 typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                 const ElementContext& elemCtx,
                 unsigned dofIdx,
                 unsigned timeIdx)
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            // TODO: let the problem do this (this is a constitutive
            // relation of which the model should be free of from the
            // abstraction POV!)
            const Evaluation& base =
                Toolbox::max(0.0001,
                             intQuants.porosity()
                             * intQuants.fluidState().saturation(phaseIdx));
            tortuosity_[phaseIdx] =
                1.0 / (intQuants.porosity() * intQuants.porosity())
                * Toolbox::pow(base, 7.0/3.0);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                diffusionCoefficient_[phaseIdx][compIdx] =
                    FluidSystem::diffusionCoefficient(fluidState,
                                                      paramCache,
                                                      phaseIdx,
                                                      compIdx);
            }
        }
    }

private:
    Evaluation tortuosity_[numPhases];
    Evaluation diffusionCoefficient_[numPhases][numComponents];
};

/*!
 * \ingroup Diffusion
 * \class Opm::DiffusionExtensiveQuantities
 *
 * \brief Provides the quantities required to calculate diffusive mass fluxes.
 */
template <class TypeTag, bool enableDiffusion>
class DiffusionExtensiveQuantities;

/*!
 * \copydoc Opm::DiffusionExtensiveQuantities
 */
template <class TypeTag>
class DiffusionExtensiveQuantities<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        the diffusive mass fluxes.
     */
    void update_(const ElementContext& elemCtx OPM_UNUSED,
                 unsigned faceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED)
    {}

    template <class Context, class FluidState>
    void updateBoundary_(const Context& context OPM_UNUSED,
                         unsigned bfIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED,
                         const FluidState& fluidState OPM_UNUSED)
    {}

public:
    /*!
     * \brief The the gradient of the mole fraction times the face normal.
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Evaluation& moleFractionGradientNormal(unsigned phaseIdx OPM_UNUSED,
                                                 unsigned compIdx OPM_UNUSED) const
    {
        throw std::logic_error("The method moleFractionGradient() does not "
                               "make sense if diffusion is disabled.");
    }

    /*!
     * \brief The effective diffusion coeffcient of a component in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Evaluation& effectiveDiffusionCoefficient(unsigned phaseIdx OPM_UNUSED,
                                                    unsigned compIdx OPM_UNUSED) const
    {
        throw std::logic_error("The method effectiveDiffusionCoefficient() "
                               "does not make sense if diffusion is disabled.");
    }
};

/*!
 * \copydoc Opm::DiffusionExtensiveQuantities
 */
template <class TypeTag>
class DiffusionExtensiveQuantities<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        the diffusive mass fluxes.
     */
    void update_(const ElementContext& elemCtx, unsigned faceIdx, unsigned timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Opm::MoleFractionCallback<TypeTag> moleFractionCallback(elemCtx);

        const auto& face = elemCtx.stencil(timeIdx).interiorFace(faceIdx);
        const auto& normal = face.normal();
        const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, timeIdx);

        const auto& intQuantsInside = elemCtx.intensiveQuantities(extQuants.interiorIndex(), timeIdx);
        const auto& intQuantsOutside = elemCtx.intensiveQuantities(extQuants.exteriorIndex(), timeIdx);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            moleFractionCallback.setPhaseIndex(phaseIdx);
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFractionCallback.setComponentIndex(compIdx);

                DimEvalVector moleFractionGradient(0.0);
                gradCalc.calculateGradient(moleFractionGradient,
                                           elemCtx,
                                           faceIdx,
                                           moleFractionCallback);

                moleFractionGradientNormal_[phaseIdx][compIdx] = 0.0;
                for (unsigned i = 0; i < normal.size(); ++i)
                    moleFractionGradientNormal_[phaseIdx][compIdx] +=
                        normal[i]*moleFractionGradient[i];
                Opm::Valgrind::CheckDefined(moleFractionGradientNormal_[phaseIdx][compIdx]);

                // use the arithmetic average for the effective
                // diffusion coefficients.
                effectiveDiffusionCoefficient_[phaseIdx][compIdx] =
                    (intQuantsInside.effectiveDiffusionCoefficient(phaseIdx, compIdx)
                     +
                     intQuantsOutside.effectiveDiffusionCoefficient(phaseIdx, compIdx))
                    / 2;

                Opm::Valgrind::CheckDefined(effectiveDiffusionCoefficient_[phaseIdx][compIdx]);
            }
        }
    }

    template <class Context, class FluidState>
    void updateBoundary_(const Context& context,
                         unsigned bfIdx,
                         unsigned timeIdx,
                         const FluidState& fluidState)
    {
        const auto& stencil = context.stencil(timeIdx);
        const auto& face = stencil.boundaryFace(bfIdx);

        const auto& elemCtx = context.elementContext();
        unsigned insideScvIdx = face.interiorIndex();
        const auto& insideScv = stencil.subControlVolume(insideScvIdx);

        const auto& intQuantsInside = elemCtx.intensiveQuantities(insideScvIdx, timeIdx);
        const auto& fluidStateInside = intQuantsInside.fluidState();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec = face.integrationPos();
        distVec -= context.element().geometry().global(insideScv.localGeometry().center());

        Scalar dist = distVec * face.normal();

        // if the following assertation triggers, the center of the
        // center of the interior SCV was not inside the element!
        assert(dist > 0);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                // calculate mole fraction gradient using two-point
                // gradients
                moleFractionGradientNormal_[phaseIdx][compIdx] =
                    (fluidState.moleFraction(phaseIdx, compIdx)
                     -
                     fluidStateInside.moleFraction(phaseIdx, compIdx))
                    / dist;
                Opm::Valgrind::CheckDefined(moleFractionGradientNormal_[phaseIdx][compIdx]);

                // use effective diffusion coefficients of the interior finite
                // volume.
                effectiveDiffusionCoefficient_[phaseIdx][compIdx] =
                    intQuantsInside.effectiveDiffusionCoefficient(phaseIdx, compIdx);

                Opm::Valgrind::CheckDefined(effectiveDiffusionCoefficient_[phaseIdx][compIdx]);
            }
        }
    }

public:
    /*!
     * \brief The the gradient of the mole fraction times the face normal.
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Evaluation& moleFractionGradientNormal(unsigned phaseIdx, unsigned compIdx) const
    { return moleFractionGradientNormal_[phaseIdx][compIdx]; }

    /*!
     * \brief The effective diffusion coeffcient of a component in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    const Evaluation& effectiveDiffusionCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return effectiveDiffusionCoefficient_[phaseIdx][compIdx]; }

private:
    Evaluation moleFractionGradientNormal_[numPhases][numComponents];
    Evaluation effectiveDiffusionCoefficient_[numPhases][numComponents];
};

} // namespace Opm

#endif
