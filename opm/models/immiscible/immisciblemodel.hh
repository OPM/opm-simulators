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
 * \copydoc Opm::ImmiscibleModel
 */
#ifndef EWOMS_IMMISCIBLE_MODEL_HH
#define EWOMS_IMMISCIBLE_MODEL_HH

#include <opm/material/densead/Math.hpp>
#include "immiscibleproperties.hh"
#include "immiscibleindices.hh"
#include "immiscibleextensivequantities.hh"
#include "immiscibleprimaryvariables.hh"
#include "immiscibleintensivequantities.hh"
#include "immiscibleratevector.hh"
#include "immiscibleboundaryratevector.hh"
#include "immisciblelocalresidual.hh"

#include <opm/models/common/multiphasebasemodel.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/SinglePhaseFluidSystem.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class ImmiscibleModel;
}

namespace Opm::Properties {

// Create new type tags
namespace TTag {
//! The generic type tag for problems using the immiscible multi-phase model
struct ImmiscibleModel { using InheritsFrom = std::tuple<VtkEnergy, MultiPhaseBaseModel>; };

//! The type tag for single-phase immiscible problems
struct ImmiscibleSinglePhaseModel { using InheritsFrom = std::tuple<ImmiscibleModel>; };

//! The type tag for two-phase immiscible problems
struct ImmiscibleTwoPhaseModel { using InheritsFrom = std::tuple<ImmiscibleModel>; };
} // end namespace TTag

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleLocalResidual<TypeTag>; };

//! the Model property
template<class TypeTag>
struct Model<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleModel<TypeTag>; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleBoundaryRateVector<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmisciblePrimaryVariables<TypeTag>; };

//! the IntensiveQuantities property
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleIntensiveQuantities<TypeTag>; };

//! the ExtensiveQuantities property
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleExtensiveQuantities<TypeTag>; };

//! The indices required by the isothermal immiscible multi-phase model
template<class TypeTag>
struct Indices<TypeTag, TTag::ImmiscibleModel> { using type = Opm::ImmiscibleIndices<TypeTag, /*PVOffset=*/0>; };

//! Disable the energy equation by default
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::ImmiscibleModel> { static constexpr bool value = false; };

/////////////////////
// set slightly different properties for the single-phase case
/////////////////////

//! The fluid system to use by default
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ImmiscibleSinglePhaseModel>
{ private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Fluid = GetPropType<TypeTag, Properties::Fluid>;
public:
    using type = Opm::SinglePhaseFluidSystem<Scalar , Fluid>;
};

template<class TypeTag>
struct Fluid<TypeTag, TTag::ImmiscibleSinglePhaseModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> >;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
template<class TypeTag>
struct VtkWriteSaturations<TypeTag, TTag::ImmiscibleSinglePhaseModel> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteMobilities<TypeTag, TTag::ImmiscibleSinglePhaseModel> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteRelativePermeabilities<TypeTag, TTag::ImmiscibleSinglePhaseModel> { static constexpr bool value = false; };

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::ImmiscibleTwoPhaseModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> >;
};

template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::ImmiscibleTwoPhaseModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> >;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ImmiscibleTwoPhaseModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;

public:
    using type = Opm::TwoPhaseImmiscibleFluidSystem<Scalar, WettingPhase, NonwettingPhase>;
};


} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup ImmiscibleModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. By default, the standard multi-phase Darcy
 * approach is used to determine the velocity, i.e.
 * \f[
 * \mathbf{v}_\alpha =
 * - \frac{k_{r\alpha}}{\mu_\alpha}
 * \mathbf{K}\left(\mathbf{grad}\, p_\alpha -
 *                 \varrho_{\alpha} \mathbf{g} \right) \;,
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
 * \frac{\partial\;\phi S_\alpha \rho_\alpha }{\partial t}
 * - \mathrm{div} \left\{ \rho_\alpha \mathbf{v}_\alpha  \right\}
 * - q_\alpha = 0 \;.
 * \f]
 *
 * The model uses the following primary variables:
 * - The pressure \f$p_0\f$ in Pascal of the phase with the lowest index
 * - The saturations \f$S_\alpha\f$ of the \f$M - 1\f$ phases that
 *   exhibit the lowest indices
 * - The absolute temperature \f$T\f$ in Kelvin if energy is conserved
 *   via the energy equation
 */
template <class TypeTag>
class ImmiscibleModel
    : public Opm::MultiPhaseBaseModel<TypeTag>
{
    using ParentType = Opm::MultiPhaseBaseModel<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { numComponents = FluidSystem::numComponents };



    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;

public:
    ImmiscibleModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "immiscible"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;

        if (pvIdx == Indices::pressure0Idx) {
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        }
        else if (Indices::saturation0Idx <= pvIdx
                 && pvIdx < Indices::saturation0Idx + numPhases - 1) {
            unsigned phaseIdx = pvIdx - Indices::saturation0Idx;
            oss << "saturation_" << FluidSystem::phaseName(phaseIdx);
        }
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

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
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
        size_t nDof = this->numTotalDof();
        for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
            if (this->isLocalDof(dofIdx)) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                break;
            }
        }
    }

    /*!
     * \copydetails FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        assert(referencePressure_ > 0);

        Scalar tmp = EnergyModule::primaryVarWeight(asImp_(), globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        if (Indices::pressure0Idx == pvIdx) {
            return 10 / referencePressure_;
        }
        return 1.0;
    }

    /*!
     * \copydetails FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(asImp_(), globalDofIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

#ifndef NDEBUG
        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(compIdx <= numPhases);
#endif

        // make all kg equal
        return 1.0;
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }

private:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    mutable Scalar referencePressure_;
};
} // namespace Opm

#endif
