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
 * \copydoc Opm::FlashModel
 */
#ifndef OPM_PTFLASH_MODEL_HH
#define OPM_PTFLASH_MODEL_HH

#include <opm/material/constraintsolvers/PTFlash.hpp>

#include <opm/material/densead/Math.hpp>

#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/models/flash/flashboundaryratevector.hh>
#include <opm/models/flash/flashextensivequantities.hh>
#include <opm/models/flash/flashproperties.hh>
#include <opm/models/flash/flashratevector.hh>

#include <opm/models/io/vtkcompositionmodule.hpp>
#include <opm/models/io/vtkdiffusionmodule.hpp>
#include <opm/models/io/vtkenergymodule.hpp>
#include <opm/models/io/vtkptflashmodule.hpp>

#include <opm/models/ptflash/flashindices.hh>
#include <opm/models/ptflash/flashintensivequantities.hh>
#include <opm/models/ptflash/flashlocalresidual.hh>
#include <opm/models/ptflash/flashnewtonmethod.hh>
#include <opm/models/ptflash/flashparameters.hh>
#include <opm/models/ptflash/flashprimaryvariables.hh>

#include <sstream>
#include <string>

namespace Opm {

template <class TypeTag>
class FlashModel;

}

namespace Opm::Properties {

namespace TTag {
//! The type tag for the isothermal single phase problems
struct FlashModel { using InheritsFrom = std::tuple<MultiPhaseBaseModel>; };
} // namespace TTag

//! Use the FlashLocalResidual function for the flash model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashLocalResidual<TypeTag>; };

//! Use the PT flash specific newton method for the flash model
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashNewtonMethod<TypeTag>; };

//! Use the Pt flash solver by default
template<class TypeTag>
struct FlashSolver<TypeTag, TTag::FlashModel>
{
    using type = Opm::PTFlash<GetPropType<TypeTag, Properties::Scalar>,
                             GetPropType<TypeTag, Properties::FluidSystem>>;
};

//! the Model property
template<class TypeTag>
struct Model<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashModel<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashPrimaryVariables<TypeTag>; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashBoundaryRateVector<TypeTag>; };

//! the IntensiveQuantities property
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashIntensiveQuantities<TypeTag>; };

//! the ExtensiveQuantities property
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashExtensiveQuantities<TypeTag>; };

//! The indices required by the flash-baseed isothermal compositional model
template<class TypeTag>
struct Indices<TypeTag, TTag::FlashModel>
{ using type = Opm::FlashIndices<TypeTag, /*PVIdx=*/0>; };

// disable molecular diffusion by default
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlashModel>
{ static constexpr bool value = false; };

//! Disable the energy equation by default
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlashModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableWater<TypeTag, TTag::MultiPhaseBaseModel> 
{ static constexpr int value = GetPropType<TypeTag, Properties::FluidSystem>::waterEnabled; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief A compositional multi-phase model based on flash-calculations
 *
 * This model assumes a flow of \f$M \geq 1\f$ fluid phases
 * \f$\alpha\f$, each of which is assumed to be a mixture \f$N \geq
 * M\f$ chemical species (denoted by the upper index \f$\kappa\f$).
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[
 * \mathbf{v}_\alpha =
 * - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 * \left(\mathbf{grad}\, p_\alpha
 *       - \varrho_{\alpha} \mathbf{g} \right) \;,
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
 * To determine the quanties that occur in the equations above, this
 * model uses <i>flash calculations</i>. A flash solver starts with
 * the total mass or molar mass per volume for each component and,
 * calculates the compositions, saturation and pressures of all
 * phases at a given temperature. For this the flash solver has to use
 * some model assumptions internally. Here a constant pressure, constant temperature,
 * two-phase flash calculation method is used.
 *
 */
template <class TypeTag>
class FlashModel
    : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };


    using EnergyModule = Opm::EnergyModule<TypeTag, enableEnergy>;

public:
    explicit FlashModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Opm::VtkCompositionModule<TypeTag>::registerParameters();
        Opm::VtkPTFlashModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Opm::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();

        Parameters::Register<Parameters::FlashTolerance<Scalar>>
            ("The maximum tolerance for the flash solver to "
             "consider the solution converged");
        Parameters::Register<Parameters::FlashVerbosity>
            ("Flash solver verbosity level");
        Parameters::Register<Parameters::FlashTwoPhaseMethod>
            ("Method for solving vapor-liquid composition. Available options include: "
             "ssi, newton, ssi+newton");

        Parameters::SetDefault<Parameters::FlashTolerance<Scalar>>(1.e-8);
        Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);

        // since thermodynamic hints are basically free if the cache for intensive quantities is
        // enabled, and this model usually shows quite a performance improvment if they are
        // enabled, let's enable them by default.
        Parameters::SetDefault<Parameters::EnableThermodynamicHints>(true);
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        const std::string& tmp = EnergyModule::primaryVarName(pvIdx);
        if (!tmp.empty())
            return tmp;

        std::ostringstream oss;
        if (Indices::z0Idx <= pvIdx && pvIdx < Indices::z0Idx + numComponents - 1)
            oss << "z_," << FluidSystem::componentName(/*compIdx=*/pvIdx - Indices::z0Idx);
        else if (pvIdx==Indices::pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(0);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(unsigned eqIdx) const
    {
        const std::string& tmp = EnergyModule::eqName(eqIdx);
        if (!tmp.empty())
            return tmp;

        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx
                                                     + numComponents) {
            unsigned compIdx = eqIdx - Indices::conti0EqIdx;
            oss << "continuity^" << FluidSystem::componentName(compIdx);
        }
        else
            assert(false);

        return oss.str();
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which are meaningful for the model
        this->addOutputModule(new Opm::VtkCompositionModule<TypeTag>(this->simulator_));
        this->addOutputModule(new Opm::VtkPTFlashModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Opm::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }
};

} // namespace Opm

#endif
