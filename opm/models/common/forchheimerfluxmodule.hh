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
 * \brief This file contains the necessary classes to calculate the
 *        volumetric fluxes out of a pressure potential gradient using the
 *        Forchhheimer approach.
 */
#ifndef EWOMS_FORCHHEIMER_FLUX_MODULE_HH
#define EWOMS_FORCHHEIMER_FLUX_MODULE_HH

#include "darcyfluxmodule.hh"

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Opm {
template <class TypeTag>
class ForchheimerIntensiveQuantities;

template <class TypeTag>
class ForchheimerExtensiveQuantities;

template <class TypeTag>
class ForchheimerBaseProblem;

/*!
 * \ingroup FluxModules
 * \brief Specifies a flux module which uses the Forchheimer relation.
 */
template <class TypeTag>
struct ForchheimerFluxModule
{
    using FluxIntensiveQuantities = ForchheimerIntensiveQuantities<TypeTag>;
    using FluxExtensiveQuantities = ForchheimerExtensiveQuantities<TypeTag>;
    using FluxBaseProblem = ForchheimerBaseProblem<TypeTag>;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    {}
};

/*!
 * \ingroup FluxModules
 * \brief Provides the defaults for the parameters required by the
 *        Forchheimer velocity approach.
 */
template <class TypeTag>
class ForchheimerBaseProblem
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    /*!
     * \brief Returns the Ergun coefficient.
     *
     * The Ergun coefficient is a measure how much the velocity is
     * reduced by turbolence. It is a quantity that does not depend on
     * the fluid phase but only on the porous medium in question. A
     * value of 0 means that the velocity is not influenced by
     * turbolence.
     */
    template <class Context>
    Scalar ergunCoefficient(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    {
        throw std::logic_error("Not implemented: Problem::ergunCoefficient()");
    }

    /*!
     * \brief Returns the ratio between the phase mobility
     *        \f$k_{r,\alpha}\f$ and its passability
     *        \f$\eta_{r,\alpha}\f$ for a given fluid phase
     *        \f$\alpha\f$.
     *
     * The passability coefficient specifies the influence of the
     * other fluid phases on the turbolent behaviour of a given fluid
     * phase. By default it is equal to the relative permeability. The
     * mobility to passability ratio is the inverse of phase' the viscosity.
     */
    template <class Context>
    Evaluation mobilityPassabilityRatio(Context& context,
                                        unsigned spaceIdx,
                                        unsigned timeIdx,
                                        unsigned phaseIdx) const
    {
        return 1.0 / context.intensiveQuantities(spaceIdx, timeIdx).fluidState().viscosity(phaseIdx);
    }
};

/*!
 * \ingroup FluxModules
 * \brief Provides the intensive quantities for the Forchheimer module
 */
template <class TypeTag>
class ForchheimerIntensiveQuantities
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

public:
    /*!
     * \brief Returns the Ergun coefficient.
     *
     * The Ergun coefficient is a measure how much the velocity is
     * reduced by turbolence. A value of 0 means that it is not
     * influenced.
     */
    const Evaluation& ergunCoefficient() const
    { return ergunCoefficient_; }

    /*!
     * \brief Returns the passability of a phase.
     */
    const Evaluation& mobilityPassabilityRatio(unsigned phaseIdx) const
    { return mobilityPassabilityRatio_[phaseIdx]; }

protected:
    void update_(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        const auto& problem = elemCtx.problem();
        ergunCoefficient_ = problem.ergunCoefficient(elemCtx, dofIdx, timeIdx);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobilityPassabilityRatio_[phaseIdx] =
                problem.mobilityPassabilityRatio(elemCtx,
                                                 dofIdx,
                                                 timeIdx,
                                                 phaseIdx);
    }

private:
    Evaluation ergunCoefficient_;
    Evaluation mobilityPassabilityRatio_[numPhases];
};

/*!
 * \ingroup FluxModules
 * \brief Provides the Forchheimer flux module
 *
 * The commonly used Darcy relation looses its validity for Reynolds numbers \f$ Re <
 * 1\f$.  If one encounters flow velocities in porous media above this threshold, the
 * Forchheimer approach can be used. Like the Darcy approach, it is a relation of with
 * the fluid velocity in terms of the gradient of pressure potential.  However, this
 * relation is not linear (as in the Darcy case) any more.
 *
 * Therefore, the Newton scheme is used to solve the Forchheimer equation. This velocity
 * is then used like the Darcy velocity e.g. by the local residual.
 *
 * Note that for Reynolds numbers above \f$\approx 500\f$ the standard Forchheimer
 * relation also looses it's validity.
 *
 * The Forchheimer equation is given by the following relation:
 *
 * \f[
  \nabla p_\alpha - \rho_\alpha \vec{g} =
  - \frac{\mu_\alpha}{k_{r,\alpha}} K^{-1}\vec{v}_\alpha
  - \frac{\rho_\alpha C_E}{\eta_{r,\alpha}} \sqrt{K}^{-1}
  \left| \vec{v}_\alpha \right| \vec{v}_\alpha
 \f]
 *
 * Where \f$C_E\f$ is the modified Ergun parameter and \f$\eta_{r,\alpha}\f$ is the
 * passability which is given by a closure relation (usually it is assumed to be
 * identical to the relative permeability). To avoid numerical problems, the relation
 * implemented by this class multiplies both sides with \f$\frac{k_{r_alpha}}{mu} K\f$,
 * so we get
 *
 * \f[
  \frac{k_{r_alpha}}{mu} K \left( \nabla p_\alpha - \rho_\alpha \vec{g}\right) =
  - \vec{v}_\alpha
  - \frac{\rho_\alpha C_E}{\eta_{r,\alpha}}  \frac{k_{r_alpha}}{mu} \sqrt{K}
  \left| \vec{v}_\alpha \right| \vec{v}_\alpha
 \f]

 */
template <class TypeTag>
class ForchheimerExtensiveQuantities
    : public DarcyExtensiveQuantities<TypeTag>
{
    using DarcyExtQuants = DarcyExtensiveQuantities<TypeTag>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using Toolbox = MathToolbox<Evaluation>;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimEvalVector = Dune::FieldVector<Evaluation, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimEvalMatrix = Dune::FieldMatrix<Evaluation, dimWorld, dimWorld>;

public:
    /*!
     * \brief Return the Ergun coefficent at the face's integration point.
     */
    const Evaluation& ergunCoefficient() const
    { return ergunCoefficient_; }

    /*!
     * \brief Return the ratio of the mobility divided by the passability at the face's
     *        integration point for a given fluid phase.
     *
     * Usually, that's the inverse of the viscosity.
     */
    Evaluation& mobilityPassabilityRatio(unsigned phaseIdx) const
    { return mobilityPassabilityRatio_[phaseIdx]; }

protected:
    void calculateGradients_(const ElementContext& elemCtx,
                             unsigned faceIdx,
                             unsigned timeIdx)
    {
        DarcyExtQuants::calculateGradients_(elemCtx, faceIdx, timeIdx);

        auto focusDofIdx = elemCtx.focusDofIndex();
        unsigned i = static_cast<unsigned>(this->interiorDofIdx_);
        unsigned j = static_cast<unsigned>(this->exteriorDofIdx_);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(j, timeIdx);

        // calculate the square root of the intrinsic permeability
        assert(isDiagonal_(this->K_));
        sqrtK_ = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            sqrtK_[dimIdx] = std::sqrt(this->K_[dimIdx][dimIdx]);

        // obtain the Ergun coefficient. Lacking better ideas, we use its the arithmetic mean.
        if (focusDofIdx == i) {
            ergunCoefficient_ =
                (intQuantsIn.ergunCoefficient() +
                 getValue(intQuantsEx.ergunCoefficient()))/2;
        }
        else if (focusDofIdx == j)
            ergunCoefficient_ =
                (getValue(intQuantsIn.ergunCoefficient()) +
                 intQuantsEx.ergunCoefficient())/2;
        else
            ergunCoefficient_ =
                (getValue(intQuantsIn.ergunCoefficient()) +
                 getValue(intQuantsEx.ergunCoefficient()))/2;

        // obtain the mobility to passability ratio for each phase.
        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            unsigned upIdx = static_cast<unsigned>(this->upstreamIndex_(phaseIdx));
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            if (focusDofIdx == upIdx) {
                density_[phaseIdx] =
                    up.fluidState().density(phaseIdx);
                mobilityPassabilityRatio_[phaseIdx] =
                    up.mobilityPassabilityRatio(phaseIdx);
            }
            else {
                density_[phaseIdx] =
                    getValue(up.fluidState().density(phaseIdx));
                mobilityPassabilityRatio_[phaseIdx] =
                    getValue(up.mobilityPassabilityRatio(phaseIdx));
            }
        }
    }

    template <class FluidState>
    void calculateBoundaryGradients_(const ElementContext& elemCtx,
                                     unsigned boundaryFaceIdx,
                                     unsigned timeIdx,
                                     const FluidState& fluidState)
    {
        DarcyExtQuants::calculateBoundaryGradients_(elemCtx,
                                                    boundaryFaceIdx,
                                                    timeIdx,
                                                    fluidState);

        auto focusDofIdx = elemCtx.focusDofIndex();
        unsigned i = static_cast<unsigned>(this->interiorDofIdx_);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);

        // obtain the Ergun coefficient. Because we are on the boundary here, we will
        // take the Ergun coefficient of the interior
        if (focusDofIdx == i)
            ergunCoefficient_ = intQuantsIn.ergunCoefficient();
        else
            ergunCoefficient_ = getValue(intQuantsIn.ergunCoefficient());

        // calculate the square root of the intrinsic permeability
        assert(isDiagonal_(this->K_));
        sqrtK_ = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            sqrtK_[dimIdx] = std::sqrt(this->K_[dimIdx][dimIdx]);

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            if (focusDofIdx == i) {
                density_[phaseIdx] = intQuantsIn.fluidState().density(phaseIdx);
                mobilityPassabilityRatio_[phaseIdx] = intQuantsIn.mobilityPassabilityRatio(phaseIdx);
            }
            else {
                density_[phaseIdx] =
                    getValue(intQuantsIn.fluidState().density(phaseIdx));
                mobilityPassabilityRatio_[phaseIdx] =
                    getValue(intQuantsIn.mobilityPassabilityRatio(phaseIdx));
            }
        }
    }

    /*!
     * \brief Calculate the volumetric fluxes of all phases
     *
     * The pressure potentials and upwind directions must already be
     * determined before calling this method!
     */
    void calculateFluxes_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        auto focusDofIdx = elemCtx.focusDofIndex();
        auto i = asImp_().interiorIndex();
        auto j = asImp_().exteriorIndex();
        const auto& intQuantsI = elemCtx.intensiveQuantities(i, timeIdx);
        const auto& intQuantsJ = elemCtx.intensiveQuantities(j, timeIdx);

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const auto& normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        // obtain the Ergun coefficient from the intensive quantity object. Until a
        // better method comes along, we use arithmetic averaging.
        if (focusDofIdx == i)
            ergunCoefficient_ =
                (intQuantsI.ergunCoefficient() +
                 getValue(intQuantsJ.ergunCoefficient())) / 2;
        else if (focusDofIdx == j)
            ergunCoefficient_ =
                (getValue(intQuantsI.ergunCoefficient()) +
                 intQuantsJ.ergunCoefficient()) / 2;
        else
            ergunCoefficient_ =
                (getValue(intQuantsI.ergunCoefficient()) +
                 getValue(intQuantsJ.ergunCoefficient())) / 2;

        ///////////////
        // calculate the weights of the upstream and the downstream control volumes
        ///////////////
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                this->filterVelocity_[phaseIdx] = 0.0;
                this->volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            calculateForchheimerFlux_(phaseIdx);

            this->volumeFlux_[phaseIdx] = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++ dimIdx)
                this->volumeFlux_[phaseIdx] +=
                    this->filterVelocity_[phaseIdx][dimIdx]*normal[dimIdx];
        }
    }

    /*!
     * \brief Calculate the volumetric flux rates of all phases at the domain boundary
     */
    void calculateBoundaryFluxes_(const ElementContext& elemCtx,
                                  unsigned bfIdx,
                                  unsigned timeIdx)
    {
        const auto& boundaryFace = elemCtx.stencil(timeIdx).boundaryFace(bfIdx);
        const auto& normal = boundaryFace.normal();

        ///////////////
        // calculate the weights of the upstream and the downstream degrees of freedom
        ///////////////
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                this->filterVelocity_[phaseIdx] = 0.0;
                this->volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            calculateForchheimerFlux_(phaseIdx);

            this->volumeFlux_[phaseIdx] = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                this->volumeFlux_[phaseIdx] +=
                    this->filterVelocity_[phaseIdx][dimIdx]*normal[dimIdx];
        }
    }

    void calculateForchheimerFlux_(unsigned phaseIdx)
    {
        // initial guess: filter velocity is zero
        DimEvalVector& velocity = this->filterVelocity_[phaseIdx];
        velocity = 0.0;

        // the change of velocity between two consecutive Newton iterations
        DimEvalVector deltaV(1e5);
        // the function value that is to be minimized of the equation that is to be
        // fulfilled
        DimEvalVector residual;
        // derivative of equation that is to be solved
        DimEvalMatrix gradResid;

        // search by means of the Newton method for a root of Forchheimer equation
        unsigned newtonIter = 0;
        while (deltaV.one_norm() > 1e-11) {
            if (newtonIter >= 50)
                throw NumericalIssue("Could not determine Forchheimer velocity within "
                                     +std::to_string(newtonIter)+" iterations");
            ++newtonIter;

            // calculate the residual and its Jacobian matrix
            gradForchheimerResid_(residual, gradResid, phaseIdx);

            // newton method
            gradResid.solve(deltaV, residual);
            velocity -= deltaV;
        }
    }

    void forchheimerResid_(DimEvalVector& residual, unsigned phaseIdx) const
    {
        const DimEvalVector& velocity = this->filterVelocity_[phaseIdx];

        // Obtaining the upstreamed quantities
        const auto& mobility = this->mobility_[phaseIdx];
        const auto& density = density_[phaseIdx];
        const auto& mobilityPassabilityRatio = mobilityPassabilityRatio_[phaseIdx];

        // optain the quantites for the integration point
        const auto& pGrad = this->potentialGrad_[phaseIdx];

        // residual = v_\alpha
        residual = velocity;

        // residual += mobility_\alpha K(\grad p_\alpha - \rho_\alpha g)
        // -> this->K_.usmv(mobility, pGrad, residual);
        assert(isDiagonal_(this->K_));
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++ dimIdx)
            residual[dimIdx] += mobility*pGrad[dimIdx]*this->K_[dimIdx][dimIdx];

        // Forchheimer turbulence correction:
        //
        // residual +=
        //   \rho_\alpha
        //   * mobility_\alpha
        //   * C_E / \eta_{r,\alpha}
        //   * abs(v_\alpha) * sqrt(K)*v_\alpha
        //
        // -> sqrtK_.usmv(density*mobilityPassabilityRatio*ergunCoefficient_*velocity.two_norm(),
        //                velocity,
        //                residual);
        Evaluation absVel = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            absVel += velocity[dimIdx]*velocity[dimIdx];
        // the derivatives of the square root of 0 are undefined, so we must guard
        // against this case
        if (absVel <= 0.0)
            absVel = 0.0;
        else
            absVel = Toolbox::sqrt(absVel);
        const auto& alpha = density*mobilityPassabilityRatio*ergunCoefficient_*absVel;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            residual[dimIdx] += sqrtK_[dimIdx]*alpha*velocity[dimIdx];
        Valgrind::CheckDefined(residual);
    }

    void gradForchheimerResid_(DimEvalVector& residual,
                               DimEvalMatrix& gradResid,
                               unsigned phaseIdx)
    {
        // TODO (?) use AD for this.
        DimEvalVector& velocity = this->filterVelocity_[phaseIdx];
        forchheimerResid_(residual, phaseIdx);

        Scalar eps = 1e-11;
        DimEvalVector tmp;
        for (unsigned i = 0; i < dimWorld; ++i) {
            Scalar coordEps = std::max(eps, Toolbox::scalarValue(velocity[i]) * (1 + eps));
            velocity[i] += coordEps;
            forchheimerResid_(tmp, phaseIdx);
            tmp -= residual;
            tmp /= coordEps;
            gradResid[i] = tmp;
            velocity[i] -= coordEps;
        }
    }

    /*!
     * \brief Check whether all off-diagonal entries of a tensor are zero.
     *
     * \param K the tensor that is to be checked.
     * \return True iff all off-diagonals are zero.
     *
     */
    bool isDiagonal_(const DimMatrix& K) const
    {
        for (unsigned i = 0; i < dimWorld; i++) {
            for (unsigned j = 0; j < dimWorld; j++) {
                if (i == j)
                    continue;

                if (std::abs(K[i][j]) > 1e-25)
                    return false;
            }
        }
        return true;
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimVector sqrtK_;

    // Ergun coefficient of all phases at the integration point
    Evaluation ergunCoefficient_;

    // Passability of all phases at the integration point
    Evaluation mobilityPassabilityRatio_[numPhases];

    // Density of all phases at the integration point
    Evaluation density_[numPhases];
};

} // namespace Opm

#endif
