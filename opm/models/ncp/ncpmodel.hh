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
 * \copydoc Opm::NcpModel
 */
#ifndef EWOMS_NCP_MODEL_HH
#define EWOMS_NCP_MODEL_HH

#include <opm/material/densead/Math.hpp>

#include "ncpproperties.hh"
#include "ncplocalresidual.hh"
#include "ncpextensivequantities.hh"
#include "ncpprimaryvariables.hh"
#include "ncpboundaryratevector.hh"
#include "ncpratevector.hh"
#include "ncpintensivequantities.hh"
#include "ncpnewtonmethod.hh"
#include "ncpindices.hh"

#include <opm/models/common/multiphasebasemodel.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace Opm {
template <class TypeTag>
class NcpModel;
}

namespace Opm::Properties {

namespace TTag {
/*!
 * \brief Define the type tag for the compositional NCP model.
 */
struct NcpModel { using InheritsFrom = std::tuple<VtkDiffusion,
                                                  VtkEnergy,
                                                  VtkComposition,
                                                  MultiPhaseBaseModel>; };
} // namespace TTag

//! Use the Ncp local jacobian operator for the compositional NCP model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NcpModel> { using type = Opm::NcpLocalResidual<TypeTag>; };

//! Use the Ncp specific newton method for the compositional NCP model
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::NcpModel> { using type = Opm::NcpNewtonMethod<TypeTag>; };

//! the Model property
template<class TypeTag>
struct Model<TypeTag, TTag::NcpModel> { using type = Opm::NcpModel<TypeTag>; };

//! The type of the base base class for actual problems
template<class TypeTag>
struct BaseProblem<TypeTag, TTag::NcpModel> { using type = Opm::MultiPhaseBaseProblem<TypeTag>; };

//! Disable the energy equation by default
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::NcpModel> { static constexpr bool value = false; };

//! disable diffusion by default
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::NcpModel> { static constexpr bool value = false; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::NcpModel> { using type = Opm::NcpRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::NcpModel> { using type = Opm::NcpBoundaryRateVector<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::NcpModel> { using type = Opm::NcpPrimaryVariables<TypeTag>; };

//! the IntensiveQuantities property
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::NcpModel> { using type = Opm::NcpIntensiveQuantities<TypeTag>; };

//! the ExtensiveQuantities property
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::NcpModel> { using type = Opm::NcpExtensiveQuantities<TypeTag>; };

//! The indices required by the compositional NCP model
template<class TypeTag>
struct Indices<TypeTag, TTag::NcpModel> { using type = Opm::NcpIndices<TypeTag, 0>; };

//! The unmodified weight for the pressure primary variable
template<class TypeTag>
struct NcpPressureBaseWeight<TypeTag, TTag::NcpModel>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
//! The weight for the saturation primary variables
template<class TypeTag>
struct NcpSaturationsBaseWeight<TypeTag, TTag::NcpModel>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
//! The unmodified weight for the fugacity primary variables
template<class TypeTag>
struct NcpFugacitiesBaseWeight<TypeTag, TTag::NcpModel>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0e-6;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup NcpModel
 *
 * \brief A compositional multi-phase model based on non-linear
 *        complementarity functions.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[
 *   \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *   \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;,
 * \f]
 * although the actual approach which is used can be specified via the
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * template<class TypeTag>
struct FluxModule<TypeTag, TTag::MyProblemTypeTag> { using type = Opm::ForchheimerFluxModule<TypeTag>; };
 * \endcode
 *
 * The core of the model is the conservation mass of each component by
 * means of the equation
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * For the missing \f$M\f$ model assumptions, the model uses
 * non-linear complementarity functions. These are based on the
 * observation that if a fluid phase is not present, the sum of the
 * mole fractions of this fluid phase is smaller than \f$1\f$, i.e.
 * \f[ \forall \alpha: S_\alpha = 0 \implies \sum_\kappa
 * x_\alpha^\kappa \leq 1 \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be non-negative:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0
 *\f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[
 * \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0
 * \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$. Such functions have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The pressure of the first phase \f$p_1\f$
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */
template <class TypeTag>
class NcpModel
    : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    using Toolbox = Opm::MathToolbox<Evaluation>;

    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;
    using DiffusionModule = Opm::DiffusionModule<TypeTag, enableDiffusion>;

public:
    NcpModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        DiffusionModule::registerParameters();
        EnergyModule::registerParameters();

        // register runtime parameters of the VTK output modules
        Opm::VtkCompositionModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Opm::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        minActivityCoeff_.resize(this->numGridDof());
        std::fill(minActivityCoeff_.begin(), minActivityCoeff_.end(), 1.0);
    }

    void adaptGrid()
    {
        ParentType::adaptGrid();
        minActivityCoeff_.resize(this->numGridDof());
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "ncp"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (saturation0Idx <= pvIdx && pvIdx < saturation0Idx + (numPhases - 1))
            oss << "saturation_" << FluidSystem::phaseName(/*phaseIdx=*/pvIdx - saturation0Idx);
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents)
            oss << "fugacity^" << FluidSystem::componentName(pvIdx - fugacity0Idx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(unsigned eqIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

        std::ostringstream oss;
        if (conti0EqIdx <= eqIdx && eqIdx < conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(eqIdx - conti0EqIdx);
        else if (ncp0EqIdx <= eqIdx && eqIdx < ncp0EqIdx + numPhases)
            oss << "ncp_" << FluidSystem::phaseName(/*phaseIdx=*/eqIdx - ncp0EqIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        // find the a reference pressure. The first degree of freedom
        // might correspond to non-interior entities which would lead
        // to an undefined value, so we have to iterate...
        for (unsigned dofIdx = 0; dofIdx < this->numGridDof(); ++ dofIdx) {
            if (this->isLocalDof(dofIdx)) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                break;
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::updatePVWeights
     */
    void updatePVWeights(const ElementContext& elemCtx) const
    {
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                minActivityCoeff_[globalIdx][compIdx] = 1e100;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const auto& fs = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState();

                    minActivityCoeff_[globalIdx][compIdx] =
                        std::min(minActivityCoeff_[globalIdx][compIdx],
                                 Toolbox::value(fs.fugacityCoefficient(phaseIdx, compIdx))
                                 * Toolbox::value(fs.pressure(phaseIdx)));
                    Opm::Valgrind::CheckDefined(minActivityCoeff_[globalIdx][compIdx]);
                }
                if (minActivityCoeff_[globalIdx][compIdx] <= 0)
                    throw Opm::NumericalIssue("The minumum activity coefficient for component "+std::to_string(compIdx)
                                                +" on DOF "+std::to_string(globalIdx)+" is negative or zero!");
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        Scalar result;
        if (tmp > 0)
            // energy related quantity
            result = tmp;
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents) {
            // component fugacity
            unsigned compIdx = pvIdx - fugacity0Idx;
            assert(compIdx <= numComponents);

            Opm::Valgrind::CheckDefined(minActivityCoeff_[globalDofIdx][compIdx]);
            static const Scalar fugacityBaseWeight =
                getPropValue<TypeTag, Properties::NcpFugacitiesBaseWeight>();
            result = fugacityBaseWeight / minActivityCoeff_[globalDofIdx][compIdx];
        }
        else if (Indices::pressure0Idx == pvIdx) {
            static const Scalar pressureBaseWeight = getPropValue<TypeTag, Properties::NcpPressureBaseWeight>();
            result = pressureBaseWeight / referencePressure_;
        }
        else {
#ifndef NDEBUG
            unsigned phaseIdx = pvIdx - saturation0Idx;
            assert(phaseIdx < numPhases - 1);
#endif

            // saturation
            static const Scalar saturationsBaseWeight =
                getPropValue<TypeTag, Properties::NcpSaturationsBaseWeight>();
            result = saturationsBaseWeight;
        }

        assert(std::isfinite(result));
        assert(result > 0);

        return result;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            // an energy related equation
            return tmp;
        // an NCP
        else if (ncp0EqIdx <= eqIdx && eqIdx < Indices::ncp0EqIdx + numPhases)
            return 1.0;

        // a mass conservation equation
        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(compIdx <= numComponents);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \brief Returns the smallest activity coefficient of a component for the
     *        most current solution at a vertex.
     *
     * \param globalDofIdx The global index of the vertex (i.e. finite volume) of interest.
     * \param compIdx The index of the component of interest.
     */
    Scalar minActivityCoeff(unsigned globalDofIdx, unsigned compIdx) const
    { return minActivityCoeff_[globalDofIdx][compIdx]; }

    /*!
     * \internal
     */
    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        this->addOutputModule(new Opm::VtkCompositionModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Opm::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }

    mutable Scalar referencePressure_;
    mutable std::vector<ComponentVector> minActivityCoeff_;
};

} // namespace Opm

#endif
