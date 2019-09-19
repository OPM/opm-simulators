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
#ifndef EWOMS_FLASH_MODEL_HH
#define EWOMS_FLASH_MODEL_HH

#include <opm/material/densead/Math.hpp>

#include "flashproperties.hh"
#include "flashprimaryvariables.hh"
#include "flashlocalresidual.hh"
#include "flashratevector.hh"
#include "flashboundaryratevector.hh"
#include "flashintensivequantities.hh"
#include "flashextensivequantities.hh"
#include "flashindices.hh"

#include <opm/models/common/multiphasebasemodel.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class FlashModel;
}

BEGIN_PROPERTIES

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(FlashModel, INHERITS_FROM(MultiPhaseBaseModel,
                                       VtkComposition,
                                       VtkEnergy,
                                       VtkDiffusion));

//! Use the FlashLocalResidual function for the flash model
SET_TYPE_PROP(FlashModel, LocalResidual,
              Opm::FlashLocalResidual<TypeTag>);

//! Use the NCP flash solver by default
SET_TYPE_PROP(FlashModel, FlashSolver,
              Opm::NcpFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
                            typename GET_PROP_TYPE(TypeTag, FluidSystem)>);

//! Let the flash solver choose its tolerance by default
SET_SCALAR_PROP(FlashModel, FlashTolerance, -1.0);

//! the Model property
SET_TYPE_PROP(FlashModel, Model, Opm::FlashModel<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(FlashModel, PrimaryVariables, Opm::FlashPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(FlashModel, RateVector, Opm::FlashRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(FlashModel, BoundaryRateVector, Opm::FlashBoundaryRateVector<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(FlashModel, IntensiveQuantities, Opm::FlashIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(FlashModel, ExtensiveQuantities, Opm::FlashExtensiveQuantities<TypeTag>);

//! The indices required by the flash-baseed isothermal compositional model
SET_TYPE_PROP(FlashModel, Indices, Opm::FlashIndices<TypeTag, /*PVIdx=*/0>);

// The updates of intensive quantities tend to be _very_ expensive for this
// model, so let's try to minimize the number of required ones
SET_BOOL_PROP(FlashModel, EnableIntensiveQuantityCache, true);

// since thermodynamic hints are basically free if the cache for intensive quantities is
// enabled, and this model usually shows quite a performance improvment if they are
// enabled, let's enable them by default.
SET_BOOL_PROP(FlashModel, EnableThermodynamicHints, true);

// disable molecular diffusion by default
SET_BOOL_PROP(FlashModel, EnableDiffusion, false);

//! Disable the energy equation by default
SET_BOOL_PROP(FlashModel, EnableEnergy, false);

END_PROPERTIES

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
 * SET_TYPE_PROP(MyProblemTypeTag, FluxModule, Opm::ForchheimerFluxModule<TypeTag>);
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
 * calculates the compositions, saturations and pressures of all
 * phases at a given temperature. For this the flash solver has to use
 * some model assumptions internally. (Often these are the same
 * primary variable switching or NCP assumptions as used by the other
 * fully implicit compositional multi-phase models provided by eWoms.)
 *
 * Using flash calculations for the flow model has some disadvantages:
 * - The accuracy of the flash solver needs to be sufficient to
 *   calculate the parital derivatives using numerical differentiation
 *   which are required for the Newton scheme.
 * - Flash calculations tend to be quite computationally expensive and
 *   are often numerically unstable.
 *
 * It is thus adviced to increase the target tolerance of the Newton
 * scheme or a to use type for scalar values which exhibits higher
 * precision than the standard \c double (e.g. \c quad) if this model
 * ought to be used.
 *
 * The model uses the following primary variables:
 * - The total molar concentration of each component:
 *   \f$c^\kappa = \sum_\alpha S_\alpha x_\alpha^\kappa \rho_{mol, \alpha}\f$
 * - The absolute temperature $T$ in Kelvins if the energy equation enabled.
 */
template <class TypeTag>
class FlashModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };


    typedef Opm::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    FlashModel(Simulator& simulator)
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

        if (enableDiffusion)
            Opm::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FlashTolerance,
                             "The maximum tolerance for the flash solver to "
                             "consider the solution converged");
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "flash"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        const std::string& tmp = EnergyModule::primaryVarName(pvIdx);
        if (tmp != "")
            return tmp;

        std::ostringstream oss;
        if (Indices::cTot0Idx <= pvIdx && pvIdx < Indices::cTot0Idx
                                                  + numComponents)
            oss << "c_tot," << FluidSystem::componentName(/*compIdx=*/pvIdx
                                                          - Indices::cTot0Idx);
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
        if (tmp != "")
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

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        if (tmp > 0)
            return tmp;

        unsigned compIdx = pvIdx - Indices::cTot0Idx;

        // make all kg equal. also, divide the weight of all total
        // compositions by 100 to make the relative errors more
        // comparable to the ones of the other models (at 10% porosity
        // the medium is fully saturated with water at atmospheric
        // conditions if 100 kg/m^3 are present!)
        return FluidSystem::molarMass(compIdx) / 100.0;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            return tmp;

        unsigned compIdx = eqIdx - Indices::conti0EqIdx;

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which are meaningful for the model
        this->addOutputModule(new Opm::VtkCompositionModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Opm::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }
};

} // namespace Opm

#endif
