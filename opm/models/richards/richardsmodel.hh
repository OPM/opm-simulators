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
 * \copydoc Opm::RichardsModel
 */
#ifndef EWOMS_RICHARDS_MODEL_HH
#define EWOMS_RICHARDS_MODEL_HH

#include <opm/material/densead/Math.hpp>

#include "richardsproperties.hh"
#include "richardsindices.hh"
#include "richardslocalresidual.hh"
#include "richardsextensivequantities.hh"
#include "richardsratevector.hh"
#include "richardsboundaryratevector.hh"
#include "richardsprimaryvariables.hh"
#include "richardsintensivequantities.hh"
#include "richardsnewtonmethod.hh"

#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class RichardsModel;
}

namespace Opm::Properties {

// Create new type tags
namespace TTag {
//! The type tag for problems discretized using the Richards model
struct Richards { using InheritsFrom = std::tuple<MultiPhaseBaseModel>; };
} // end namespace TTag

//! By default, assume that the first phase is the liquid one
template<class TypeTag>
struct LiquidPhaseIndex<TypeTag, TTag::Richards> { static constexpr int value = 0; };

//! By default, assume that the non-liquid phase is gaseos
template<class TypeTag>
struct GasPhaseIndex<TypeTag, TTag::Richards> { static constexpr int value = 1 - getPropValue<TypeTag, Properties::LiquidPhaseIndex>(); };

/*!
 * \brief By default, assume that component which the liquid is made of has
 *        the same index as the liquid phase.
 *
 * This is a convention which works for most fluid systems shipped
 * with eWoms by default, but it cannot generally correct because the
 * liquid can be composed of different components. (e.g., do you
 * prefer Ethanol of H2O??)
 */
template<class TypeTag>
struct LiquidComponentIndex<TypeTag, TTag::Richards> { static constexpr int value = getPropValue<TypeTag, Properties::LiquidPhaseIndex>(); };

//! By default, assume that the gas component is the other than the liquid one
template<class TypeTag>
struct GasComponentIndex<TypeTag, TTag::Richards> { static constexpr int value = 1 - getPropValue<TypeTag, Properties::LiquidComponentIndex>(); };

//! The local residual operator
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Richards> { using type = Opm::RichardsLocalResidual<TypeTag>; };

//! The global model used
template<class TypeTag>
struct Model<TypeTag, TTag::Richards> { using type = Opm::RichardsModel<TypeTag>; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::Richards> { using type = Opm::RichardsRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::Richards> { using type = Opm::RichardsBoundaryRateVector<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::Richards> { using type = Opm::RichardsPrimaryVariables<TypeTag>; };

//! The class for the intensive quantities
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::Richards> { using type = Opm::RichardsIntensiveQuantities<TypeTag>; };

//! The class for the quantities required for the flux calculation
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::Richards> { using type = Opm::RichardsExtensiveQuantities<TypeTag>; };

//! The class of the Newton method
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::Richards> { using type = Opm::RichardsNewtonMethod<TypeTag>; };

//! The class with all index definitions for the model
template<class TypeTag>
struct Indices<TypeTag, TTag::Richards> { using type = Opm::RichardsIndices; };

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the Richards model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
template<class TypeTag>
struct WettingFluid<TypeTag, TTag::Richards>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> >;
};

/*!
 * \brief The non-wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the Richards model to work because the
 * Richards model does not conserve the non-wetting phase.
 */
template<class TypeTag>
struct NonWettingFluid<TypeTag, TTag::Richards>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::GasPhase<Scalar, Opm::NullComponent<Scalar> >;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingFluid and NonWettingFluid properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Richards>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingFluid = GetPropType<TypeTag, Properties::WettingFluid>;
    using NonWettingFluid = GetPropType<TypeTag, Properties::NonWettingFluid>;

public:
    using type = Opm::TwoPhaseImmiscibleFluidSystem<Scalar, WettingFluid, NonWettingFluid>;
};


} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup RichardsModel
 *
 * \brief This model implements a variant of the Richards equation for
 *        quasi-twophase flow.
 *
 * In the unsaturated zone, Richards' equation is frequently used to
 * approximate the water distribution above the groundwater level. It
 * can be derived from the two-phase equations, i.e.
 * \f[
 * \frac{\partial\;\phi S_\alpha \rho_\alpha}{\partial t}
 * -
 * \mathrm{div} \left\{
 * \rho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}\; \mathbf{K}\;
 * \mathbf{grad}\left[
 * p_\alpha - g\rho_\alpha
 * \right]
 * \right\}
 * =
 * q_\alpha,
 * \f]
 * where \f$\alpha \in \{w, n\}\f$ is the index of the fluid phase,
 * \f$\rho_\alpha\f$ is the fluid density, \f$S_\alpha\f$ is the fluid
 * saturation, \f$\phi\f$ is the porosity of the soil,
 * \f$k_{r\alpha}\f$ is the relative permeability for the fluid,
 * \f$\mu_\alpha\f$ is the fluid's dynamic viscosity, \f$\mathbf{K}\f$
 * is the intrinsic permeability tensor, \f$p_\alpha\f$ is the fluid
 * phase pressure and \f$g\f$ is the potential of the gravity field.
 *
 * In contrast to the "full" two-phase model, the Richards model
 * assumes that the non-wetting fluid is gas and that it thus exhibits
 * a much lower viscosity than the (liquid) wetting phase. (This
 * assumption is quite realistic in many applications: For example, at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. Using this reasoning, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitely large compared to the same
 * term of the liquid phase. This implies that the pressure of the gas
 * phase is equivalent to the static pressure distribution and that
 * therefore, mass conservation only needs to be considered for the
 * liquid phase.
 *
 * The model thus choses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 * \f[
 * S_w = p_c^{-1}(p_n - p_w)
 * \f]
 * holds, where \f$p_n\f$ is a reference pressure given by the
 * problem's \c referencePressure() method. Nota bene, that the last
 * step assumes that the capillary pressure-saturation curve can be
 * uniquely inverted, i.e. it is not possible to set the capillary
 * pressure to zero if the Richards model ought to be used!
 */
template <class TypeTag>
class RichardsModel
    : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

     static const unsigned numPhases = FluidSystem::numPhases;
     static const unsigned numComponents = FluidSystem::numComponents;

    static const unsigned liquidPhaseIdx = getPropValue<TypeTag, Properties::LiquidPhaseIndex>();
    static const unsigned gasPhaseIdx = getPropValue<TypeTag, Properties::GasPhaseIndex>();

    static const unsigned liquidCompIdx = getPropValue<TypeTag, Properties::LiquidComponentIndex>();
    static const unsigned gasCompIdx = getPropValue<TypeTag, Properties::GasComponentIndex>();


    // some consistency checks
    static_assert(numPhases == 2,
                  "Exactly two fluids are required for this model");
    static_assert(numComponents == 2,
                  "Exactly two components are required for this model");
    static_assert(liquidPhaseIdx != gasPhaseIdx,
                  "The liquid and the gas phases must be different");
    static_assert(liquidCompIdx != gasCompIdx,
                  "The liquid and the gas components must be different");

public:
    RichardsModel(Simulator& simulator)
        : ParentType(simulator)
    {
        // the liquid phase must be liquid, the gas phase must be
        // gaseous. Think about it!
        assert(FluidSystem::isLiquid(liquidPhaseIdx));
        assert(!FluidSystem::isLiquid(gasPhaseIdx));
    }

    /*!
     * \copydoc FvBaseDiscretization::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "richards"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::ostringstream oss;
        if (pvIdx == Indices::pressureWIdx)
            oss << "pressure_" << FluidSystem::phaseName(liquidPhaseIdx);
        else
            assert(0);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(unsigned eqIdx) const
    {
        std::ostringstream oss;
        if (eqIdx == Indices::contiEqIdx)
            oss << "continuity_" << FluidSystem::phaseName(liquidPhaseIdx);
        else
            assert(0);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned, unsigned pvIdx) const
    {
        if (Indices::pressureWIdx == pvIdx) {
            return 10 / referencePressure_;
        }

        return 1;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned, [[maybe_unused]] unsigned eqIdx) const
    {
        assert((eqIdx - Indices::contiEqIdx) <= FluidSystem::numPhases);

        // make all kg equal
        return 1.0;
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
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressureWIdx];
                break;
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::phaseIsConsidered
     */
    bool phaseIsConsidered(unsigned phaseIdx) const
    { return phaseIdx == liquidPhaseIdx; }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();
    }

    mutable Scalar referencePressure_;
};
} // namespace Opm

#endif
