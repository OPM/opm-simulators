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
 * \copydoc Opm::ObstacleProblem
 */
#ifndef EWOMS_OBSTACLE_PROBLEM_HH
#define EWOMS_OBSTACLE_PROBLEM_HH

#include <opm/models/ncp/ncpproperties.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <iostream>

namespace Opm {
template <class TypeTag>
class ObstacleProblem;
}

namespace Opm::Properties {

namespace TTag {
struct ObstacleBaseProblem {};
}

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ObstacleBaseProblem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ObstacleBaseProblem> { using type = Opm::ObstacleProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ObstacleBaseProblem>
{ using type = Opm::H2ON2FluidSystem<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::ObstacleBaseProblem>
{
private:
    // define the material law
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialTraits = Opm::TwoPhaseMaterialTraits<Scalar,
                                                       /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                                       /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    using EffMaterialLaw = Opm::LinearMaterial<MaterialTraits>;

public:
    using type = Opm::EffToAbsLaw<EffMaterialLaw>;
};

// Set the thermal conduction law
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::ObstacleBaseProblem>
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
struct SolidEnergyLaw<TypeTag, TTag::ObstacleBaseProblem>
{ using type = Opm::ConstantSolidHeatCapLaw<GetPropType<TypeTag, Properties::Scalar>>; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::ObstacleBaseProblem> { static constexpr bool value = true; };

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::ObstacleBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e4;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::ObstacleBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 250;
};

// The default DGF file to load
template<class TypeTag>
struct GridFile<TypeTag, TTag::ObstacleBaseProblem> { static constexpr auto value = "./data/obstacle_24x16.dgf"; };

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where liquid water is first stopped by a
 *        low-permeability lens and then seeps though it.
 *
 * Liquid water is injected by using of a free-flow condition on the
 * lower right of the domain. This water level then raises until
 * hydrostatic pressure is reached. On the left of the domain, a
 * rectangular obstacle with \f$10^3\f$ lower permeability than the
 * rest of the domain first stops the for a while until it seeps
 * through it.
 *
 * The domain is sized 60m times 40m and consists of two media, a
 * moderately permeable soil (\f$ K_0=10e-12 m^2\f$) and an obstacle
 * at \f$[10; 20]m \times [0; 35]m \f$ with a lower permeablility of
 * \f$ K_1=K_0/1000\f$.
 *
 * Initially the whole domain is filled by nitrogen, the temperature
 * is \f$20^\circ C\f$ for the whole domain. The gas pressure is
 * initially 1 bar, at the inlet of the liquid water on the right side
 * it is 2 bar.
 *
 * The boundary is no-flow except on the lower 10 meters of the left
 * and the right boundary where a free flow condition is assumed.
 */
template <class TypeTag>
class ObstacleProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using ThermalConductionLawParams = GetPropType<TypeTag, Properties::ThermalConductionLawParams>;
    using SolidEnergyLawParams = GetPropType<TypeTag, Properties::SolidEnergyLawParams>;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>(),
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dimWorld>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ObstacleProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        unsigned nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        unsigned np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

        // intrinsic permeabilities
        coarseK_ = this->toDimMatrix_(1e-12);
        fineK_ = this->toDimMatrix_(1e-15);

        // the porosity
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.0);
        fineMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the linear law, i.e. minimum and maximum
        // pressures
        fineMaterialParams_.setPcMinSat(liquidPhaseIdx, 0.0);
        fineMaterialParams_.setPcMaxSat(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setPcMinSat(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setPcMaxSat(liquidPhaseIdx, 0.0);

        /*
        // entry pressures for Brooks-Corey
        fineMaterialParams_.setEntryPressure(5e3);
        coarseMaterialParams_.setEntryPressure(1e3);

        // Brooks-Corey shape parameters
        fineMaterialParams_.setLambda(2);
        coarseMaterialParams_.setLambda(2);
        */

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law of thermal conduction
        computeThermalCondParams_(fineThermalCondParams_, finePorosity_);
        computeThermalCondParams_(coarseThermalCondParams_, coarsePorosity_);

        // assume constant volumetric heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();

        initFluidStates_();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

        // Calculate storage terms of the individual phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            this->model().globalPhaseStorage(phaseStorage, phaseIdx);

            if (this->gridView().comm().rank() == 0) {
                std::cout << "Storage in " << FluidSystem::phaseName(phaseIdx)
                          << "Phase: [" << phaseStorage << "]"
                          << "\n"  << std::flush;
            }
        }

        // Calculate total storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage total: [" << storage << "]"
                      << "\n"  << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "obstacle"
            << "_" << Model::name();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     *
     * This problem simply assumes a constant temperature.
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix&
    intrinsicPermeability(const Context& context,
                          unsigned spaceIdx,
                          unsigned timeIdx) const
    {
        if (isFineMaterial_(context.pos(spaceIdx, timeIdx)))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context,
                    unsigned spaceIdx,
                    unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& context,
                      unsigned spaceIdx,
                      unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

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
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionParams(const Context& context,
                         unsigned spaceIdx,
                         unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineThermalCondParams_;
        return coarseThermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);

        if (onInlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, inletFluidState_);
        else if (onOutlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, outletFluidState_);
        else
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context,
                 unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(outletFluidState_, matParams);
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
    { rate = 0.0; }

    //! \}

private:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     */
    bool isFineMaterial_(const GlobalPosition& pos) const
    { return 10 <= pos[0] && pos[0] <= 20 && 0 <= pos[1] && pos[1] <= 35; }

    bool onInlet_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10;
    }

    bool onOutlet_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10;
    }

    void initFluidStates_()
    {
        initFluidState_(inletFluidState_, coarseMaterialParams_,
                        /*isInlet=*/true);
        initFluidState_(outletFluidState_, coarseMaterialParams_,
                        /*isInlet=*/false);
    }

    template <class FluidState>
    void initFluidState_(FluidState& fs, const MaterialLawParams& matParams, bool isInlet)
    {
        unsigned refPhaseIdx;
        unsigned otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(temperature_);

        if (isInlet) {
            // only liquid on inlet
            refPhaseIdx = liquidPhaseIdx;
            otherPhaseIdx = gasPhaseIdx;

            // set liquid saturation
            fs.setSaturation(liquidPhaseIdx, 1.0);

            // set pressure of the liquid phase
            fs.setPressure(liquidPhaseIdx, 2e5);

            // set the liquid composition to pure water
            fs.setMoleFraction(liquidPhaseIdx, N2Idx, 0.0);
            fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);
        }
        else {
            // elsewhere, only gas
            refPhaseIdx = gasPhaseIdx;
            otherPhaseIdx = liquidPhaseIdx;

            // set gas saturation
            fs.setSaturation(gasPhaseIdx, 1.0);

            // set pressure of the gas phase
            fs.setPressure(gasPhaseIdx, 1e5);

            // set the gas composition to 99% nitrogen and 1% steam
            fs.setMoleFraction(gasPhaseIdx, N2Idx, 0.99);
            fs.setMoleFraction(gasPhaseIdx, H2OIdx, 0.01);
        }

        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        fs.setPressure(otherPhaseIdx, fs.pressure(refPhaseIdx)
                                      + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

        // make the fluid state consistent with local thermodynamic
        // equilibrium
        using ComputeFromReferencePhase = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        ComputeFromReferencePhase::solve(fs, paramCache, refPhaseIdx,
                                         /*setViscosity=*/true,
                                         /*setEnthalpy=*/false);
    }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    DimMatrix coarseK_;
    DimMatrix fineK_;

    Scalar coarsePorosity_;
    Scalar finePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams fineThermalCondParams_;
    ThermalConductionLawParams coarseThermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> inletFluidState_;
    Opm::CompositionalFluidState<Scalar, FluidSystem> outletFluidState_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Opm

#endif
