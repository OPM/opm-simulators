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
 * \brief Contains the classes required to consider energy as a
 *        conservation quantity in a multi-phase module.
 */
#ifndef EWOMS_ENERGY_MODULE_HH
#define EWOMS_ENERGY_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

BEGIN_PROPERTIES

NEW_PROP_TAG(Indices);
NEW_PROP_TAG(EnableEnergy);
NEW_PROP_TAG(ThermalConductionLaw);
NEW_PROP_TAG(ThermalConductionLawParams);
NEW_PROP_TAG(SolidEnergyLaw);
NEW_PROP_TAG(SolidEnergyLawParams);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup Energy
 * \brief Provides the auxiliary methods required for consideration of
 *        the energy equation.
 */
template <class TypeTag, bool enableEnergy>
class EnergyModule;

/*!
 * \copydoc Opm::EnergyModule
 */
template <class TypeTag>
class EnergyModule<TypeTag, /*enableEnergy=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Evaluation, numEq> EvalEqVector;

public:
    /*!
     * \brief Register all run-time parameters for the energy module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Returns the name of a primary variable or an empty
     *        string if the specified primary variable index does not belong to
     *        the energy module.
     */
    static std::string primaryVarName(unsigned pvIdx OPM_UNUSED)
    { return ""; }

    /*!
     * \brief Returns the name of an equation or an empty
     *        string if the specified equation index does not belong to
     *        the energy module.
     */
    static std::string eqName(unsigned eqIdx OPM_UNUSED)
    { return ""; }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     */
    static Scalar primaryVarWeight(const Model& model OPM_UNUSED,
                                   unsigned globalDofIdx OPM_UNUSED,
                                   unsigned pvIdx OPM_UNUSED)
    { return -1; }

    /*!
     * \brief Returns the relative weight of a equation of the residual.
     */
    static Scalar eqWeight(const Model& model OPM_UNUSED,
                           unsigned globalDofIdx OPM_UNUSED,
                           unsigned eqIdx OPM_UNUSED)
    { return -1; }

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables& priVars OPM_UNUSED,
                                      const FluidState& fs OPM_UNUSED)
    {}

    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class RateVector, class FluidState>
    static void setEnthalpyRate(RateVector& rateVec OPM_UNUSED,
                                const FluidState& fluidState OPM_UNUSED,
                                unsigned phaseIdx OPM_UNUSED,
                                const Evaluation& volume OPM_UNUSED)
    {}

    /*!
     * \brief Add the rate of the enthalpy flux to a rate vector.
     */
    static void setEnthalpyRate(RateVector& rateVec OPM_UNUSED,
                                const Evaluation& rate OPM_UNUSED)
    {}

    /*!
     * \brief Add the rate of the enthalpy flux to a rate vector.
     */
    static void addToEnthalpyRate(RateVector& rateVec OPM_UNUSED,
                                  const Evaluation& rate OPM_UNUSED)
    {}

    /*!
     * \brief Add the rate of the conductive energy flux to a rate vector.
     */
    static Scalar thermalConductionRate(const ExtensiveQuantities& extQuants OPM_UNUSED)
    { return 0.0; }

    /*!
     * \brief Add the energy storage term for a fluid phase to an equation
     * vector
     */
    template <class LhsEval>
    static void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage OPM_UNUSED,
                                const IntensiveQuantities& intQuants OPM_UNUSED,
                                unsigned phaseIdx OPM_UNUSED)
    {}

    /*!
     * \brief Add the energy storage term for a fluid phase to an equation
     * vector
     */
    template <class LhsEval, class Scv>
    static void addFracturePhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage OPM_UNUSED,
                                        const IntensiveQuantities& intQuants OPM_UNUSED,
                                        const Scv& scv OPM_UNUSED,
                                        unsigned phaseIdx OPM_UNUSED)
    {}

    /*!
     * \brief Add the energy storage term for the fracture part a fluid phase to an
     *        equation vector
     */
    template <class LhsEval>
    static void addSolidEnergyStorage(Dune::FieldVector<LhsEval, numEq>& storage OPM_UNUSED,
                                    const IntensiveQuantities& intQuants OPM_UNUSED)
    {}

    /*!
     * \brief Evaluates the advective energy fluxes over a face of a
     *        subcontrol volume and adds the result in the flux vector.
     *
     * This method is called by compute flux (base class)
     */
    template <class Context>
    static void addAdvectiveFlux(RateVector& flux OPM_UNUSED,
                                 const Context& context OPM_UNUSED,
                                 unsigned spaceIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED)
    {}

    /*!
     * \brief Evaluates the advective energy fluxes over a fracture
     *        which should be attributed to a face of a subcontrol
     *        volume and adds the result in the flux vector.
     */
    template <class Context>
    static void handleFractureFlux(RateVector& flux OPM_UNUSED,
                                   const Context& context OPM_UNUSED,
                                   unsigned spaceIdx OPM_UNUSED,
                                   unsigned timeIdx OPM_UNUSED)
    {}

    /*!
     * \brief Adds the diffusive energy flux to the flux vector over the face of a
     *        sub-control volume.
     *
     * This method is called by compute flux (base class)
     */
    template <class Context>
    static void addDiffusiveFlux(RateVector& flux OPM_UNUSED,
                                 const Context& context OPM_UNUSED,
                                 unsigned spaceIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED)
    {}
};

/*!
 * \copydoc Opm::EnergyModule
 */
template <class TypeTag>
class EnergyModule<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };
    enum { energyEqIdx = Indices::energyEqIdx };
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef Dune::FieldVector<Evaluation, numEq> EvalEqVector;
    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \brief Register all run-time parameters for the energy module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Returns the name of a primary variable or an empty
     *        string if the specified primary variable index does not belong to
     *        the energy module.
     */
    static std::string primaryVarName(unsigned pvIdx)
    {
        if (pvIdx == temperatureIdx)
            return "temperature";
        return "";
    }

    /*!
     * \brief Returns the name of an equation or an empty
     *        string if the specified equation index does not belong to
     *        the energy module.
     */
    static std::string eqName(unsigned eqIdx)
    {
        if (eqIdx == energyEqIdx)
            return "energy";
        return "";
    }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     */
    static Scalar primaryVarWeight(const Model& model, unsigned globalDofIdx, unsigned pvIdx)
    {
        if (pvIdx != temperatureIdx)
            return -1;

        // make the weight of the temperature primary variable inversly proportional to its value
        return std::max(1.0/1000, 1.0/model.solution(/*timeIdx=*/0)[globalDofIdx][temperatureIdx]);
    }

    /*!
     * \brief Returns the relative weight of a equation.
     */
    static Scalar eqWeight(const Model& model OPM_UNUSED,
                           unsigned globalDofIdx OPM_UNUSED,
                           unsigned eqIdx)
    {
        if (eqIdx != energyEqIdx)
            return -1;

        // approximate change of internal energy of 1kg of liquid water for a temperature
        // change of 30K
        return 1.0 / (4.184e3 * 30.0);
    }

    /*!
     * \brief Set the rate of energy flux of a rate vector.
     */
    static void setEnthalpyRate(RateVector& rateVec, const Evaluation& rate)
    { rateVec[energyEqIdx] = rate; }

    /*!
     * \brief Add the rate of the enthalpy flux to a rate vector.
     */
    static void addToEnthalpyRate(RateVector& rateVec, const Evaluation& rate)
    { rateVec[energyEqIdx] += rate; }

    /*!
     * \brief Returns the conductive energy flux for a given flux integration point.
     */
    static Evaluation thermalConductionRate(const ExtensiveQuantities& extQuants)
    { return -extQuants.temperatureGradNormal() * extQuants.thermalConductivity(); }

    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class RateVector, class FluidState>
    static void setEnthalpyRate(RateVector& rateVec,
                                const FluidState& fluidState,
                                unsigned phaseIdx,
                                const Evaluation& volume)
    {
        rateVec[energyEqIdx] =
            volume
            * fluidState.density(phaseIdx)
            * fluidState.enthalpy(phaseIdx);
    }

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables& priVars,
                                      const FluidState& fs)
    {
        priVars[temperatureIdx] = Toolbox::value(fs.temperature(/*phaseIdx=*/0));
#ifndef NDEBUG
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(std::abs(Toolbox::value(fs.temperature(/*phaseIdx=*/0))
                            - Toolbox::value(fs.temperature(phaseIdx))) < 1e-30);
        }
#endif
    }

    /*!
     * \brief Add the energy storage term for a fluid phase to an equation
     * vector
     */
    template <class LhsEval>
    static void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                                const IntensiveQuantities& intQuants,
                                unsigned phaseIdx)
    {
        const auto& fs = intQuants.fluidState();
        storage[energyEqIdx] +=
            Toolbox::template decay<LhsEval>(fs.density(phaseIdx))
            * Toolbox::template decay<LhsEval>(fs.internalEnergy(phaseIdx))
            * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
            * Toolbox::template decay<LhsEval>(intQuants.porosity());
    }

    /*!
     * \brief Add the energy storage term for a fluid phase to an equation
     * vector
     */
    template <class Scv, class LhsEval>
    static void addFracturePhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                                        const IntensiveQuantities& intQuants,
                                        const Scv& scv,
                                        unsigned phaseIdx)
    {
        const auto& fs = intQuants.fractureFluidState();
        storage[energyEqIdx] +=
            Toolbox::template decay<LhsEval>(fs.density(phaseIdx))
            * Toolbox::template decay<LhsEval>(fs.internalEnergy(phaseIdx))
            * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
            * Toolbox::template decay<LhsEval>(intQuants.fracturePorosity())
            * Toolbox::template decay<LhsEval>(intQuants.fractureVolume())/scv.volume();
    }

    /*!
     * \brief Add the energy storage term for a fluid phase to an equation
     *        vector
     */
    template <class LhsEval>
    static void addSolidEnergyStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                                    const IntensiveQuantities& intQuants)
    { storage[energyEqIdx] += Opm::decay<LhsEval>(intQuants.solidInternalEnergy()); }

    /*!
     * \brief Evaluates the advective energy fluxes for a flux integration point and adds
     *        the result in the flux vector.
     *
     * This method is called by compute flux (base class)
     */
    template <class Context>
    static void addAdvectiveFlux(RateVector& flux,
                                 const Context& context,
                                 unsigned spaceIdx,
                                 unsigned timeIdx)
    {
        const auto& extQuants = context.extensiveQuantities(spaceIdx, timeIdx);

        // advective energy flux in all phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!context.model().phaseIsConsidered(phaseIdx))
                continue;

            // intensive quantities of the upstream and the downstream DOFs
            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = context.intensiveQuantities(upIdx, timeIdx);

            flux[energyEqIdx] +=
                extQuants.volumeFlux(phaseIdx)
                * up.fluidState().enthalpy(phaseIdx)
                * up.fluidState().density(phaseIdx);
        }
    }

    /*!
     * \brief Evaluates the advective energy fluxes over a fracture which should be
     *        attributed to a face of a subcontrol volume and adds the result in the flux
     *        vector.
     */
    template <class Context>
    static void handleFractureFlux(RateVector& flux,
                                   const Context& context,
                                   unsigned spaceIdx,
                                   unsigned timeIdx)
    {
        const auto& scvf = context.stencil(timeIdx).interiorFace(spaceIdx);
        const auto& extQuants = context.extensiveQuantities(spaceIdx, timeIdx);

        // reduce the energy flux in the matrix by the half the width occupied by the
        // fracture
        flux[energyEqIdx] *=
            1 - extQuants.fractureWidth()/(2*scvf.area());

        // advective energy flux in all phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!context.model().phaseIsConsidered(phaseIdx))
                continue;

            // intensive quantities of the upstream and the downstream DOFs
            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = context.intensiveQuantities(upIdx, timeIdx);

            flux[energyEqIdx] +=
                extQuants.fractureVolumeFlux(phaseIdx)
                * up.fluidState().enthalpy(phaseIdx)
                * up.fluidState().density(phaseIdx);
        }
    }

    /*!
     * \brief Adds the diffusive energy flux to the flux vector over the face of a
     *        sub-control volume.
     *
     * This method is called by compute flux (base class)
     */
    template <class Context>
    static void addDiffusiveFlux(RateVector& flux,
                                 const Context& context,
                                 unsigned spaceIdx,
                                 unsigned timeIdx)
    {
        const auto& extQuants = context.extensiveQuantities(spaceIdx, timeIdx);

        // diffusive energy flux
        flux[energyEqIdx] +=
            - extQuants.temperatureGradNormal()
            * extQuants.thermalConductivity();
    }
};

/*!
 * \ingroup Energy
 * \class Opm::EnergyIndices
 *
 * \brief Provides the indices required for the energy equation.
 */
template <unsigned PVOffset, bool enableEnergy>
struct EnergyIndices;

/*!
 * \copydoc Opm::EnergyIndices
 */
template <unsigned PVOffset>
struct EnergyIndices<PVOffset, /*enableEnergy=*/false>
{
    //! The index of the primary variable representing temperature
    enum { temperatureIdx = -1000 };

    //! The index of the equation representing the conservation of energy
    enum { energyEqIdx = -1000 };

protected:
    enum { numEq_ = 0 };
};

/*!
 * \copydoc Opm::EnergyIndices
 */
template <unsigned PVOffset>
struct EnergyIndices<PVOffset, /*enableEnergy=*/true>
{
    //! The index of the primary variable representing temperature
    enum { temperatureIdx = PVOffset };

    //! The index of the equation representing the conservation of energy
    enum { energyEqIdx = PVOffset };

protected:
    enum { numEq_ = 1 };
};

/*!
 * \ingroup Energy
 * \class Opm::EnergyIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the energy equation.
 */
template <class TypeTag, bool enableEnergy>
class EnergyIntensiveQuantities;

/*!
 * \copydoc Opm::EnergyIntensiveQuantities
 */
template <class TypeTag>
class EnergyIntensiveQuantities<TypeTag, /*enableEnergy=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \brief Returns the volumetric internal energy \f$\mathrm{[J/(m^3]}\f$ of the
     *        solid matrix in the sub-control volume.
     */
    Evaluation solidInternalEnergy() const
    {
        throw std::logic_error("solidInternalEnergy() does not make sense for isothermal models");
    }

    /*!
     * \brief Returns the total thermal conductivity \f$\mathrm{[W/m^2 / (K/m)]}\f$ of
     *        the solid matrix in the sub-control volume.
     */
    Evaluation thermalConductivity() const
    {
        throw std::logic_error("thermalConductivity() does not make sense for isothermal models");
    }

protected:
    /*!
     * \brief Update the temperatures of the fluids of a fluid state.
     */
    template <class FluidState, class Context>
    static void updateTemperatures_(FluidState& fluidState,
                                    const Context& context,
                                    unsigned spaceIdx,
                                    unsigned timeIdx)
    {
        Scalar T = context.problem().temperature(context, spaceIdx, timeIdx);
        fluidState.setTemperature(Toolbox::createConstant(T));
    }

    /*!
     * \brief Update the quantities required to calculate
     *        energy fluxes.
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
 * \copydoc Opm::EnergyIntensiveQuantities
 */
template <class TypeTag>
class EnergyIntensiveQuantities<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLaw) ThermalConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLaw) SolidEnergyLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases = FluidSystem::numPhases };
    enum { energyEqIdx = Indices::energyEqIdx };
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef Opm::MathToolbox<Evaluation> Toolbox;

protected:
    /*!
     * \brief Update the temperatures of the fluids of a fluid state.
     */
    template <class FluidState, class Context>
    static void updateTemperatures_(FluidState& fluidState,
                                    const Context& context,
                                    unsigned spaceIdx,
                                    unsigned timeIdx)
    {
        const auto& priVars = context.primaryVars(spaceIdx, timeIdx);
        Evaluation val;
        if (std::is_same<Evaluation, Scalar>::value) // finite differences
            val = Toolbox::createConstant(priVars[temperatureIdx]);
        else {
            // automatic differentiation
            if (timeIdx == 0)
                val = Toolbox::createVariable(priVars[temperatureIdx], temperatureIdx);
            else
                val = Toolbox::createConstant(priVars[temperatureIdx]);
        }
        fluidState.setTemperature(val);
    }

    /*!
     * \brief Update the quantities required to calculate
     *        energy fluxes.
     */
    template <class FluidState>
    void update_(FluidState& fs,
                 typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                 const ElementContext& elemCtx,
                 unsigned dofIdx,
                 unsigned timeIdx)
    {
        // set the specific enthalpies of the fluids
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            fs.setEnthalpy(phaseIdx,
                           FluidSystem::enthalpy(fs, paramCache, phaseIdx));
        }

        // compute and set the volumetric internal energy of the solid phase
        const auto& problem = elemCtx.problem();
        const auto& solidEnergyParams = problem.solidEnergyLawParams(elemCtx, dofIdx, timeIdx);
        const auto& thermalCondParams = problem.thermalConductionLawParams(elemCtx, dofIdx, timeIdx);

        solidInternalEnergy_ = SolidEnergyLaw::solidInternalEnergy(solidEnergyParams, fs);
        thermalConductivity_ = ThermalConductionLaw::thermalConductivity(thermalCondParams, fs);

        Opm::Valgrind::CheckDefined(solidInternalEnergy_);
        Opm::Valgrind::CheckDefined(thermalConductivity_);
    }

public:
    /*!
     * \brief Returns the volumetric internal energy \f$\mathrm{[J/m^3]}\f$ of the
     *        solid matrix in the sub-control volume.
     */
    const Evaluation& solidInternalEnergy() const
    { return solidInternalEnergy_; }

    /*!
     * \brief Returns the total conductivity capacity \f$\mathrm{[W/m^2 / (K/m)]}\f$ of
     *        the solid matrix in the sub-control volume.
     */
    const Evaluation& thermalConductivity() const
    { return thermalConductivity_; }

private:
    Evaluation solidInternalEnergy_;
    Evaluation thermalConductivity_;
};

/*!
 * \ingroup Energy
 * \class Opm::EnergyExtensiveQuantities
 *
 * \brief Provides the quantities required to calculate energy fluxes.
 */
template <class TypeTag, bool enableEnergy>
class EnergyExtensiveQuantities;

/*!
 * \copydoc Opm::EnergyExtensiveQuantities
 */
template <class TypeTag>
class EnergyExtensiveQuantities<TypeTag, /*enableEnergy=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        energy fluxes.
     */
    void update_(const ElementContext& elemCtx OPM_UNUSED,
                 unsigned faceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED)
    {}

    template <class Context, class FluidState>
    void updateBoundary_(const Context& context OPM_UNUSED,
                         unsigned bfIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED,
                         const FluidState& fs OPM_UNUSED)
    {}

public:
    /*!
     * \brief The temperature gradient times the face normal [K m^2 / m]
     */
    Scalar temperatureGradNormal() const
    {
        throw std::logic_error("Calling temperatureGradNormal() does not make sense "
                               "for isothermal models");
    }

    /*!
     * \brief The total thermal conductivity at the face \f$\mathrm{[W/m^2 / (K/m)]}\f$
     */
    Scalar thermalConductivity() const
    {
        throw std::logic_error("Calling thermalConductivity() does not make sense for "
                               "isothermal models");
    }
};

/*!
 * \copydoc Opm::EnergyExtensiveQuantities
 */
template <class TypeTag>
class EnergyExtensiveQuantities<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        energy fluxes.
     */
    void update_(const ElementContext& elemCtx, unsigned faceIdx, unsigned timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Opm::TemperatureCallback<TypeTag> temperatureCallback(elemCtx);

        EvalDimVector temperatureGrad;
        gradCalc.calculateGradient(temperatureGrad,
                                   elemCtx,
                                   faceIdx,
                                   temperatureCallback);

        // scalar product of temperature gradient and scvf normal
        const auto& face = elemCtx.stencil(/*timeIdx=*/0).interiorFace(faceIdx);

        temperatureGradNormal_ = 0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            temperatureGradNormal_ += (face.normal()[dimIdx]*temperatureGrad[dimIdx]);

        const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, timeIdx);
        const auto& intQuantsInside = elemCtx.intensiveQuantities(extQuants.interiorIndex(), timeIdx);
        const auto& intQuantsOutside = elemCtx.intensiveQuantities(extQuants.exteriorIndex(), timeIdx);

        // arithmetic mean
        thermalConductivity_ =
            0.5 * (intQuantsInside.thermalConductivity() + intQuantsOutside.thermalConductivity());
        Opm::Valgrind::CheckDefined(thermalConductivity_);
    }

    template <class Context, class FluidState>
    void updateBoundary_(const Context& context, unsigned bfIdx, unsigned timeIdx, const FluidState& fs)
    {
        const auto& stencil = context.stencil(timeIdx);
        const auto& face = stencil.boundaryFace(bfIdx);

        const auto& elemCtx = context.elementContext();
        unsigned insideScvIdx = face.interiorIndex();
        const auto& insideScv = elemCtx.stencil(timeIdx).subControlVolume(insideScvIdx);

        const auto& intQuantsInside = elemCtx.intensiveQuantities(insideScvIdx, timeIdx);
        const auto& fsInside = intQuantsInside.fluidState();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec = face.integrationPos();
        distVec -= insideScv.geometry().center();

        Scalar tmp = 0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            tmp += distVec[dimIdx] * face.normal()[dimIdx];
        Scalar dist = tmp;

        // if the following assertation triggers, the center of the
        // center of the interior SCV was not inside the element!
        assert(dist > 0);

        // calculate the temperature gradient using two-point gradient
        // appoximation
        temperatureGradNormal_ =
            (fs.temperature(/*phaseIdx=*/0) - fsInside.temperature(/*phaseIdx=*/0)) / dist;

        // take the value for thermal conductivity from the interior finite volume
        thermalConductivity_ = intQuantsInside.thermalConductivity();
    }

public:
    /*!
     * \brief The temperature gradient times the face normal [K m^2 / m]
     */
    const Evaluation& temperatureGradNormal() const
    { return temperatureGradNormal_; }

    /*!
     * \brief The total thermal conductivity at the face \f$\mathrm{[W/m^2 /
     * (K/m)]}\f$
     */
    const Evaluation& thermalConductivity() const
    { return thermalConductivity_; }

private:
    Evaluation temperatureGradNormal_;
    Evaluation thermalConductivity_;
};

} // namespace Opm

#endif
