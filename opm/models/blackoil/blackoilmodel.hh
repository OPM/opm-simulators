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
 * \copydoc Opm::BlackOilModel
 */
#ifndef OPM_BLACK_OIL_MODEL_HPP
#define OPM_BLACK_OIL_MODEL_HPP

#include <opm/material/densead/Math.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilboundaryratevector.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>
#include <opm/models/blackoil/blackoildarcyfluxmodule.hh>
#include <opm/models/blackoil/blackoildiffusionmodule.hh>
#include <opm/models/blackoil/blackoildispersionmodule.hh>
#include <opm/models/blackoil/blackoilextbomodules.hh>
#include <opm/models/blackoil/blackoilextensivequantities.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/blackoil/blackoillocalresidual.hh>
#include <opm/models/blackoil/blackoilmicpmodules.hh>
#include <opm/models/blackoil/blackoilnewtonmethod.hpp>
#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilprimaryvariables.hh>
#include <opm/models/blackoil/blackoilproblem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/blackoil/blackoilratevector.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/models/io/vtkblackoilmodule.hpp>
#include <opm/models/io/vtkcompositionmodule.hpp>
#include <opm/models/io/vtkdiffusionmodule.hpp>

#include <cassert>
#include <istream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace Opm {

template <class TypeTag>
class BlackOilModel;

}

namespace Opm::Properties {

namespace TTag {

//! The type tag for the black-oil problems
struct BlackOilModel
{ using InheritsFrom = std::tuple<VtkBlackOilPolymer, MultiPhaseBaseModel>; };
} // namespace TTag

//! Set the local residual function
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilLocalResidual<TypeTag>; };

//! Use the black-oil specific newton method
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilNewtonMethod<TypeTag>; };

//! The Model property
template<class TypeTag>
struct Model<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilModel<TypeTag>; };

//! The Problem property
template<class TypeTag>
struct BaseProblem<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilProblem<TypeTag>; };

//! the RateVector property
template<class TypeTag>
struct RateVector<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilRateVector<TypeTag>; };

//! the BoundaryRateVector property
template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilBoundaryRateVector<TypeTag>; };

//! the PrimaryVariables property
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilPrimaryVariables<TypeTag>; };

//! the IntensiveQuantities property
template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilIntensiveQuantities<TypeTag>; };

//! the ExtensiveQuantities property
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilExtensiveQuantities<TypeTag>; };

//! Use the the velocity module which is aware of the black-oil specific model extensions
//! (i.e., the polymer and solvent extensions)
template<class TypeTag>
struct FluxModule<TypeTag, TTag::BlackOilModel>
{ using type = BlackOilDarcyFluxModule<TypeTag>; };

//! The indices required by the model
template<class TypeTag>
struct Indices<TypeTag, TTag::BlackOilModel>
{
    using type = BlackOilVariableAndEquationIndices<getPropValue<TypeTag, Properties::EnableSolvent>(),
                                                    getPropValue<TypeTag, Properties::EnableExtbo>(),
                                                    getPropValue<TypeTag, Properties::EnablePolymer>(),
                                                    getPropValue<TypeTag, Properties::EnableEnergy>(),
                                                    getPropValue<TypeTag, Properties::EnableFoam>(),
                                                    getPropValue<TypeTag, Properties::EnableBrine>(),
                                                    /*PVOffset=*/0,
                                                    getPropValue<TypeTag, Properties::EnableMICP>()>;
};

//! Set the fluid system to the black-oil fluid system by default
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BlackOilModel>
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using type = BlackOilFluidSystem<Scalar>;
};

// by default, all ECL extension modules are disabled
template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableExtbo<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnablePolymerMW<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableFoam<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableBrine<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableVapwat<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableDisgasInWater<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableSaltPrecipitation<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableMICP<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

//! By default, the blackoil model is isothermal and does not conserve energy
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

//! disable diffusion by default
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

//! disable disperison by default
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

//! by default, scale the energy equation by the inverse of the energy required to heat
//! up one kg of water by 30 Kelvin. If we conserve surface volumes, this must be divided
//! by the weight of one cubic meter of water. This is required to make the "dumb" linear
//! solvers that do not weight the components of the solutions do the right thing.
//! by default, don't scale the energy equation, i.e. assume that a reasonable linear
//! solver is used. (Not scaling it makes debugging quite a bit easier.)
template<class TypeTag>
struct BlackOilEnergyScalingFactor<TypeTag, TTag::BlackOilModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr Scalar alpha =
        getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>() ? 1000.0 : 1.0;

public:
    using type = Scalar;
    static constexpr Scalar value = 1.0/(30.0*4184.0*alpha);
};

//! similarly to the energy equation, the urea equation in MICP is scaled by a factor of 0.1
template<class TypeTag>
struct BlackOilUreaScalingFactor<TypeTag, TTag::BlackOilModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr Scalar alpha =
        getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>() ? 1000.0 : 1.0;

public:
    using type = Scalar;
    static constexpr Scalar value = 1.0/(10.0*alpha);
};

// by default, ebos formulates the conservation equations in terms of mass not surface
// volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::BlackOilModel>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilModel
 * \brief A fully-implicit black-oil flow model.
 *
 * The black-oil model is a three-phase, three-component model widely
 * used for oil reservoir simulation.  The phases are denoted by lower
 * index \f$\alpha \in \{ w, g, o \}\f$ ("water", "gas" and "oil") and
 * the components by upper index \f$\kappa \in \{ W, G, O \}\f$
 * ("Water", "Gas" and "Oil"). The model assumes partial miscibility:
 *
 * - Water and the gas phases are immisicible and are assumed to be
 *   only composed of the water and gas components respectively-
 * - The oil phase is assumed to be a mixture of the gas and the oil
 *  components.
 *
 * The densities of the phases are determined by so-called
 * <i>formation volume factors</i>:
 *
 * \f[
 * B_\alpha := \frac{\varrho_\alpha(1\,\text{bar})}{\varrho_\alpha(p_\alpha)}
 * \f]
 *
 * Since the gas and water phases are assumed to be immiscible, this
 * is sufficient to calculate their density. For the formation volume
 * factor of the the oil phase \f$B_o\f$ determines the density of
 * *saturated* oil, i.e. the density of the oil phase if some gas
 * phase is present.
 *
 * The composition of the oil phase is given by the <i>gas dissolution factor</i>
 * \f$R_s\f$, which defined as the volume of gas at atmospheric pressure that is
 * dissolved in a given amount of oil at reservoir pressure:
 *
 * \f[
 * R_s := \frac{\varrho_{o}^G}{\varrho_o^O}\;.
 * \f]
 *
 * This allows to calculate all quantities required for the
 * mass-conservation equations for each component, i.e.
 *
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;,
 * \f]
 * where \f$\mathrm{v}_\alpha\f$ is the filter velocity of the phase
 * \f$\alpha\f$.
 *
 * By default \f$\mathrm{v}_\alpha\f$ is determined by using the
 * standard multi-phase Darcy approach, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *\left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * template<class TypeTag>
struct FluxModule<TypeTag, TTag::MyProblemTypeTag> { using type = Opm::ForchheimerFluxModule<TypeTag>; };
 * \endcode
 *
 * The primary variables used by this model are:
 * - The pressure of the phase with the lowest index
 * - The two saturations of the phases with the lowest indices
 */
template<class TypeTag >
class BlackOilModel
    : public MultiPhaseBaseModel<TypeTag>
{
public:
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

private:
    using Implementation = GetPropType<TypeTag, Properties::Model>;
    using ParentType = MultiPhaseBaseModel<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = FluidSystem::numComponents };
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;
    static constexpr bool waterEnabled = Indices::waterEnabled;

    using SolventModule = BlackOilSolventModule<TypeTag>;
    using ExtboModule = BlackOilExtboModule<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using EnergyModule = BlackOilEnergyModule<TypeTag>;
    using DiffusionModule = BlackOilDiffusionModule<TypeTag, enableDiffusion>;
    using DispersionModule = BlackOilDispersionModule<TypeTag, enableDispersion>;
    using MICPModule = BlackOilMICPModule<TypeTag>;

public:
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;

    explicit BlackOilModel(Simulator& simulator)
        : ParentType(simulator)
    {
        eqWeights_.resize(numEq, 1.0);
    }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        SolventModule::registerParameters();
        ExtboModule::registerParameters();
        PolymerModule::registerParameters();
        EnergyModule::registerParameters();
        DiffusionModule::registerParameters();
        MICPModule::registerParameters();

        // register runtime parameters of the VTK output modules
        VtkBlackOilModule<TypeTag>::registerParameters();
        VtkCompositionModule<TypeTag>::registerParameters();
        VtkDiffusionModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "blackoil"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        if (pvIdx == Indices::waterSwitchIdx) {
            return "water_switching";
        }
        else if (pvIdx == Indices::pressureSwitchIdx) {
            return "pressure_switching";
        }
        else if (static_cast<int>(pvIdx) == Indices::compositionSwitchIdx) {
            return "composition_switching";
        }
        else if (SolventModule::primaryVarApplies(pvIdx)) {
            return SolventModule::primaryVarName(pvIdx);
        }
        else if (ExtboModule::primaryVarApplies(pvIdx)) {
            return ExtboModule::primaryVarName(pvIdx);
        }
        else if (PolymerModule::primaryVarApplies(pvIdx)) {
            return PolymerModule::primaryVarName(pvIdx);
        }
        else if (EnergyModule::primaryVarApplies(pvIdx)) {
            return EnergyModule::primaryVarName(pvIdx);
        }
        else {
            throw std::logic_error("Invalid primary variable index");
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(int eqIdx) const
    {
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents) {
            std::ostringstream oss;
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
            return oss.str();
        }
        else if (SolventModule::eqApplies(eqIdx)) {
            return SolventModule::eqName(eqIdx);
        }
        else if (ExtboModule::eqApplies(eqIdx)) {
            return ExtboModule::eqName(eqIdx);
        }
        else if (PolymerModule::eqApplies(eqIdx)) {
            return PolymerModule::eqName(eqIdx);
        }
        else if (EnergyModule::eqApplies(eqIdx)) {
            return EnergyModule::eqName(eqIdx);
        }
        else {
            throw std::logic_error("Invalid equation index");
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= this->numGridDof()) {
            return 1.0;
        }

        // saturations are always in the range [0, 1]!
        if (int(Indices::waterSwitchIdx) == int(pvIdx)) {
            return 1.0;
        }

        // oil pressures usually are in the range of 100 to 500 bars for typical oil
        // reservoirs (which is the only relevant application for the black-oil model).
        else if (int(Indices::pressureSwitchIdx) == int(pvIdx)) {
            return 1.0 / 300e5;
        }

        // deal with primary variables stemming from the solvent module
        else if (SolventModule::primaryVarApplies(pvIdx)) {
            return SolventModule::primaryVarWeight(pvIdx);
        }

        // deal with primary variables stemming from the extBO module
        else if (ExtboModule::primaryVarApplies(pvIdx)) {
            return ExtboModule::primaryVarWeight(pvIdx);
        }

        // deal with primary variables stemming from the polymer module
        else if (PolymerModule::primaryVarApplies(pvIdx)) {
            return PolymerModule::primaryVarWeight(pvIdx);
        }

        // deal with primary variables stemming from the energy module
        else if (EnergyModule::primaryVarApplies(pvIdx)) {
            return EnergyModule::primaryVarWeight(pvIdx);
        }

        // if the primary variable is either the gas saturation, Rs or Rv
        assert(int(Indices::compositionSwitchIdx) == int(pvIdx));

        switch (this->solution(0)[globalDofIdx].primaryVarsMeaningGas()) {
        case PrimaryVariables::GasMeaning::Sg: return 1.0; // gas saturation
        case PrimaryVariables::GasMeaning::Rs: return 1.0 / 250.; // gas dissolution factor
        case PrimaryVariables::GasMeaning::Rv: return 1.0 / 0.025; // oil vaporization factor
        default: throw std::logic_error("Invalid primary variable meaning flag for gas");
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= this->numGridDof()) {
            return 1.0;
        }

        return eqWeights_[eqIdx];
    }

    void setEqWeight(unsigned eqIdx, Scalar value)
    { eqWeights_[eqIdx] = value; }

    /*!
     * \brief Write the current solution for a degree of freedom to a
     *        restart file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param dof The Dune entity which's data should be serialized
     */
    template <class DofEntity>
    void serializeEntity(std::ostream& outstream, const DofEntity& dof)
    {
        const unsigned dofIdx = static_cast<unsigned>(asImp_().dofMapper().index(dof));

        // write phase state
        if (!outstream.good()) {
            throw std::runtime_error("Could not serialize degree of freedom " + std::to_string(dofIdx));
        }

        // write the primary variables
        const auto& priVars = this->solution(/*timeIdx=*/0)[dofIdx];
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << priVars[eqIdx] << " ";
        }

        // write the pseudo primary variables
        outstream << static_cast<int>(priVars.primaryVarsMeaningGas()) << " ";
        outstream << static_cast<int>(priVars.primaryVarsMeaningWater()) << " ";
        outstream << static_cast<int>(priVars.primaryVarsMeaningPressure()) << " ";

        outstream << priVars.pvtRegionIndex() << " ";

        SolventModule::serializeEntity(asImp_(), outstream, dof);
        ExtboModule::serializeEntity(asImp_(), outstream, dof);
        PolymerModule::serializeEntity(asImp_(), outstream, dof);
        EnergyModule::serializeEntity(asImp_(), outstream, dof);
    }

    /*!
     * \brief Reads the current solution variables for a degree of
     *        freedom from a restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param dof The Dune entity which's data should be deserialized
     */
    template <class DofEntity>
    void deserializeEntity(std::istream& instream,
                           const DofEntity& dof)
    {
        const unsigned dofIdx = static_cast<unsigned>(asImp_().dofMapper().index(dof));

        // read in the "real" primary variables of the DOF
        auto& priVars = this->solution(/*timeIdx=*/0)[dofIdx];
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good()) {
                throw std::runtime_error("Could not deserialize degree of freedom " + std::to_string(dofIdx));
            }
            instream >> priVars[eqIdx];
        }

        // read the pseudo primary variables
        unsigned primaryVarsMeaningGas;
        instream >> primaryVarsMeaningGas;

        unsigned primaryVarsMeaningWater;
        instream >> primaryVarsMeaningWater;

        unsigned primaryVarsMeaningPressure;
        instream >> primaryVarsMeaningPressure;

        unsigned pvtRegionIdx;
        instream >> pvtRegionIdx;

        if (!instream.good()) {
            throw std::runtime_error("Could not deserialize degree of freedom " + std::to_string(dofIdx));
        }

        SolventModule::deserializeEntity(asImp_(), instream, dof);
        ExtboModule::deserializeEntity(asImp_(), instream, dof);
        PolymerModule::deserializeEntity(asImp_(), instream, dof);
        EnergyModule::deserializeEntity(asImp_(), instream, dof);

        using PVM_G = typename PrimaryVariables::GasMeaning;
        using PVM_W = typename PrimaryVariables::WaterMeaning;
        using PVM_P = typename PrimaryVariables::PressureMeaning;
        priVars.setPrimaryVarsMeaningGas(static_cast<PVM_G>(primaryVarsMeaningGas));
        priVars.setPrimaryVarsMeaningWater(static_cast<PVM_W>(primaryVarsMeaningWater));
        priVars.setPrimaryVarsMeaningPressure(static_cast<PVM_P>(primaryVarsMeaningPressure));

        priVars.setPvtRegionIndex(pvtRegionIdx);
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        ParentType::deserialize(res);

        // set the PVT indices of the primary variables. This is also done by writing
        // them into the restart file and re-reading them, but it is better to calculate
        // them from scratch because the input could have been changed in this regard...
        ElementContext elemCtx(this->simulator_);
        for (const auto& elem : elements(this->gridView())) {
            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timIdx=*/0); ++dofIdx) {
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timIdx=*/0);
                updatePvtRegionIndex_(this->solution(/*timeIdx=*/0)[globalDofIdx],
                                      elemCtx,
                                      dofIdx,
                                      /*timeIdx=*/0);
            }
        }

        this->solution(/*timeIdx=*/1) = this->solution(/*timeIdx=*/0);
    }

/*
    // hack: this interferes with the static polymorphism trick
protected:
    friend ParentType;
    friend Discretization;
*/

    template <class Context>
    void supplementInitialSolution_(PrimaryVariables& priVars,
                                    const Context& context,
                                    unsigned dofIdx,
                                    unsigned timeIdx)
    { updatePvtRegionIndex_(priVars, context, dofIdx, timeIdx); }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which make sense for the blackoil model
        SolventModule::registerOutputModules(asImp_(), this->simulator_);
        PolymerModule::registerOutputModules(asImp_(), this->simulator_);
        EnergyModule::registerOutputModules(asImp_(), this->simulator_);
        MICPModule::registerOutputModules(asImp_(), this->simulator_);

        this->addOutputModule(std::make_unique<VtkBlackOilModule<TypeTag>>(this->simulator_));
        this->addOutputModule(std::make_unique<VtkCompositionModule<TypeTag>>(this->simulator_));

        if constexpr (enableDiffusion) {
            this->addOutputModule(std::make_unique<VtkDiffusionModule<TypeTag>>(this->simulator_));
        }
    }

private:
    std::vector<Scalar> eqWeights_;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template <class Context>
    void updatePvtRegionIndex_(PrimaryVariables& priVars,
                               const Context& context,
                               unsigned dofIdx,
                               unsigned timeIdx)
    {
        const unsigned regionIdx = context.problem().pvtRegionIndex(context, dofIdx, timeIdx);
        priVars.setPvtRegionIndex(regionIdx);
    }
};

} // namespace Opm

#endif // OPM_BLACK_OIL_MODEL_HPP
