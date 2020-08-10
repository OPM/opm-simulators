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
 * \copydoc Opm::FractureProblem
 */
#ifndef EWOMS_FRACTURE_PROBLEM_HH
#define EWOMS_FRACTURE_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
// avoid reordering of macro elements, otherwise this problem won't work
#define DISABLE_ALUGRID_SFC_ORDERING 1
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#else
#error "dune-alugrid not found!"
#endif

#include <opm/models/discretefracture/discretefracturemodel.hh>
#include <opm/models/io/dgfvanguard.hh>

#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Dnapl.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <iostream>
#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class FractureProblem;
}

namespace Opm::Properties {

// Create a type tag for the problem
// Create new type tags
namespace TTag {
struct FractureProblem { using InheritsFrom = std::tuple<DiscreteFractureModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FractureProblem>
{ using type = Dune::ALUGrid</*dim=*/2, /*dimWorld=*/2, Dune::simplex, Dune::nonconforming>; };

// Set the Vanguard property
template<class TypeTag>
struct Vanguard<TypeTag, TTag::FractureProblem> { using type = Opm::DgfVanguard<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FractureProblem> { using type = Opm::FractureProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::FractureProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::FractureProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::DNAPL<Scalar> >;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FractureProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffectiveLaw = Opm::RegularizedBrooksCorey<Traits>;
    // using EffectiveLaw = RegularizedVanGenuchten<Traits>;
    // using EffectiveLaw = LinearMaterial<Traits>;
public:
    using type = Opm::EffToAbsLaw<EffectiveLaw>;
};

// Enable the energy equation
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FractureProblem> { static constexpr bool value = true; };

// Set the thermal conduction law
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::FractureProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::SomertonThermalConductionLaw<FluidSystem, Scalar>;
};

// set the energy storage law for the solid phase
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FractureProblem>
{ using type = Opm::ConstantSolidHeatCapLaw<GetPropType<TypeTag, Properties::Scalar>>; };

// Disable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::FractureProblem> { static constexpr bool value = false; };

// For this problem, we use constraints to specify the left boundary
template<class TypeTag>
struct EnableConstraints<TypeTag, TTag::FractureProblem> { static constexpr bool value = true; };

// Set the default value for the file name of the grid
template<class TypeTag>
struct GridFile<TypeTag, TTag::FractureProblem> { static constexpr auto value = "data/fracture.art.dgf"; };

// Set the default value for the end time
template<class TypeTag>
struct EndTime<TypeTag, TTag::FractureProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3e3;
};

// Set the default value for the initial time step size
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::FractureProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;
};

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Two-phase problem which involves fractures
 *
 * The domain is initially completely saturated by the oil phase,
 * except for the left side, which is fully water saturated. Since the
 * capillary pressure in the fractures is lower than in the rock
 * matrix and the material is hydrophilic, water infiltrates through
 * the fractures and gradually pushes the oil out on the right side,
 * where the pressure is kept constant.
 */
template <class TypeTag>
class FractureProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using ThermalConductionLawParams = GetPropType<TypeTag, Properties::ThermalConductionLawParams>;
    using SolidEnergyLawParams = GetPropType<TypeTag, Properties::SolidEnergyLawParams>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    enum {
        // phase indices
        wettingPhaseIdx = MaterialLaw::wettingPhaseIdx,
        nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx,

        // number of phases
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using FluidState = Opm::ImmiscibleFluidState<Scalar, FluidSystem>;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    template <int dim>
    struct FaceLayout
    {
        bool contains(Dune::GeometryType gt)
        { return gt.dim() == dim - 1; }
    };
    using FaceMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    using FractureMapper = Opm::FractureMapper<TypeTag>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    FractureProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C

        matrixMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        matrixMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);
        fractureMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        fractureMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

#if 0 // linear
        matrixMaterialParams_.setEntryPC(0.0);
        matrixMaterialParams_.setMaxPC(2000.0);
        fractureMaterialParams_.setEntryPC(0.0);
        fractureMaterialParams_.setMaxPC(1000.0);
#endif

#if 1 // Brooks-Corey
        matrixMaterialParams_.setEntryPressure(2000);
        matrixMaterialParams_.setLambda(2.0);
        matrixMaterialParams_.setPcLowSw(1e-1);
        fractureMaterialParams_.setEntryPressure(1000);
        fractureMaterialParams_.setLambda(2.0);
        fractureMaterialParams_.setPcLowSw(5e-2);
#endif

#if 0 // van Genuchten
        matrixMaterialParams_.setVgAlpha(0.0037);
        matrixMaterialParams_.setVgN(4.7);
        fractureMaterialParams_.setVgAlpha(0.0025);
        fractureMaterialParams_.setVgN(4.7);
#endif

        matrixMaterialParams_.finalize();
        fractureMaterialParams_.finalize();

        matrixK_ = this->toDimMatrix_(1e-15);         // m^2
        fractureK_ = this->toDimMatrix_(1e5 * 1e-15); // m^2

        matrixPorosity_ = 0.10;
        fracturePorosity_ = 0.25;
        fractureWidth_ = 1e-3; // [m]

        // initialize the energy-related parameters
        initEnergyParams_(thermalConductionParams_, matrixPorosity_);
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "fracture_" << Model::name();
        return oss.str();
    }

    /*!
     * \brief Called directly after the time integration.
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        // checkConservativeness() does not include the effect of constraints, so we
        // disable it for this problem...
        //this->model().checkConservativeness();

        // Calculate storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    // \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    { return matrixK_; }

    /*!
     * \brief Intrinsic permeability of fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const DimMatrix& fractureIntrinsicPermeability(const Context& context OPM_UNUSED,
                                                   unsigned spaceIdx OPM_UNUSED,
                                                   unsigned timeIdx OPM_UNUSED) const
    { return fractureK_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return matrixPorosity_; }

    /*!
     * \brief The porosity inside the fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    Scalar fracturePorosity(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    { return fracturePorosity_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context OPM_UNUSED,
                                               unsigned spaceIdx OPM_UNUSED,
                                               unsigned timeIdx OPM_UNUSED) const
    { return matrixMaterialParams_; }

    /*!
     * \brief The parameters for the material law inside the fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const MaterialLawParams& fractureMaterialLawParams(const Context& context OPM_UNUSED,
                                                       unsigned spaceIdx OPM_UNUSED,
                                                       unsigned timeIdx OPM_UNUSED) const
    { return fractureMaterialParams_; }

    /*!
     * \brief Returns the object representating the fracture topology.
     */
    const FractureMapper& fractureMapper() const
    { return this->simulator().vanguard().fractureMapper(); }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \todo This method should get one face index instead of two
     *       vertex indices. This probably requires a new context
     *       class, though.
     *
     * \param context The execution context.
     * \param spaceIdx1 The local index of the edge's first edge.
     * \param spaceIdx2 The local index of the edge's second edge.
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar fractureWidth(const Context& context OPM_UNUSED,
                         unsigned spaceIdx1 OPM_UNUSED,
                         unsigned spaceIdx2 OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return fractureWidth_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams&
    thermalConductionLawParams(const Context& context OPM_UNUSED,
                               unsigned spaceIdx OPM_UNUSED,
                               unsigned timeIdx OPM_UNUSED) const
    { return thermalConductionParams_; }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyParams_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(pos)) {
            // on the right boundary, we impose a free-flow
            // (i.e. Dirichlet) condition
            FluidState fluidState;
            fluidState.setTemperature(temperature_);

            fluidState.setSaturation(wettingPhaseIdx, 0.0);
            fluidState.setSaturation(nonWettingPhaseIdx,
                                     1.0 - fluidState.saturation(wettingPhaseIdx));

            fluidState.setPressure(wettingPhaseIdx, 1e5);
            fluidState.setPressure(nonWettingPhaseIdx, fluidState.pressure(wettingPhaseIdx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fluidState);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fluidState.setDensity(phaseIdx,
                                      FluidSystem::density(fluidState, paramCache, phaseIdx));
                fluidState.setViscosity(phaseIdx,
                                        FluidSystem::viscosity(fluidState, paramCache, phaseIdx));
            }

            // set a free flow (i.e. Dirichlet) boundary
            values.setFreeFlow(context, spaceIdx, timeIdx, fluidState);
        }
        else
            // for the upper, lower and left boundaries, use a no-flow
            // condition (i.e. a Neumann 0 condition)
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::constraints
     */
    template <class Context>
    void constraints(Constraints& constraints, const Context& context,
                     unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (!onLeftBoundary_(pos))
            // only impose constraints adjacent to the left boundary
            return;

        unsigned globalIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        if (!fractureMapper().isFractureVertex(globalIdx)) {
            // do not impose constraints if the finite volume does
            // not contain fractures.
            return;
        }

        // if the current finite volume is on the left boundary
        // and features a fracture, specify the fracture fluid
        // state.
        FluidState fractureFluidState;
        fractureFluidState.setTemperature(temperature_ + 10.0);

        fractureFluidState.setSaturation(wettingPhaseIdx, 1.0);
        fractureFluidState.setSaturation(nonWettingPhaseIdx,
                                         1.0 - fractureFluidState.saturation(
                                                   wettingPhaseIdx));

        Scalar pCFracture[numPhases];
        MaterialLaw::capillaryPressures(pCFracture, fractureMaterialParams_,
                                        fractureFluidState);

        fractureFluidState.setPressure(wettingPhaseIdx, /*pressure=*/1.0e5);
        fractureFluidState.setPressure(nonWettingPhaseIdx,
                                       fractureFluidState.pressure(wettingPhaseIdx)
                                       + (pCFracture[nonWettingPhaseIdx]
                                          - pCFracture[wettingPhaseIdx]));

        constraints.setActive(true);
        constraints.assignNaiveFromFracture(fractureFluidState,
                                            matrixMaterialParams_);
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        FluidState fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wettingPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(nonWettingPhaseIdx, fluidState.pressure(wettingPhaseIdx));

        fluidState.setSaturation(wettingPhaseIdx, 0.0);
        fluidState.setSaturation(nonWettingPhaseIdx,
                                 1.0 - fluidState.saturation(wettingPhaseIdx));

        values.assignNaive(fluidState);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    void initEnergyParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        // assume the volumetric heat capacity of granite
        solidEnergyParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                * 2700.0); // density of granite [kg/m^3]
        solidEnergyParams_.finalize();

        Scalar lambdaGranite = 2.8; // [W / (K m)]

        // create a Fluid state which has all phases present
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(293.15);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, 1.0135e5);
        }

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar lambdaSaturated;
            if (FluidSystem::isLiquid(phaseIdx)) {
                Scalar lambdaFluid = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
                lambdaSaturated =
                    std::pow(lambdaGranite, (1 - poro))
                    + std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1 - poro));

            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
        }

        Scalar lambdaVac = std::pow(lambdaGranite, (1 - poro));
        params.setVacuumLambda(lambdaVac);
    }

    DimMatrix matrixK_;
    DimMatrix fractureK_;

    Scalar matrixPorosity_;
    Scalar fracturePorosity_;

    Scalar fractureWidth_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams matrixMaterialParams_;

    ThermalConductionLawParams thermalConductionParams_;
    SolidEnergyLawParams solidEnergyParams_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Opm

#endif // EWOMS_FRACTURE_PROBLEM_HH
