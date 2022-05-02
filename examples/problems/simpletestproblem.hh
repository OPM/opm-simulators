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
 * \copydoc Opm::simpletestproblem
 */
#ifndef EWOMS_SIMPLETEST_PROBLEM_HH
#define EWOMS_SIMPLETEST_PROBLEM_HH

#include <opm/common/Exceptions.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/BrooksCorey.hpp>
#include <opm/material/constraintsolvers/PTFlash.hpp> 
#include <opm/material/fluidsystems/ThreeComponentFluidSystem.hh>
#include <opm/material/common/Valgrind.hpp>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/ptflash/flashmodel.hh>
#include <opm/models/io/structuredgridvanguard.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>
#include <opm/simulators/linalg/parallelistlbackend.hh>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Opm {
template <class TypeTag>
class SimpleTest;
} // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct SimpleTest {};
} // end namespace TTag

// declare the "simpletest" problem specify property tags
template <class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct SimulationName { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct EpisodeLength { using type = UndefinedProperty;};

template <class TypeTag, class MyTypeTag>
struct Initialpressure { using type = UndefinedProperty;};

// Set the grid type: --->1D
template <class TypeTag>
struct Grid<TypeTag, TTag::SimpleTest> { using type = Dune::YaspGrid</*dim=*/2>; };

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::SimpleTest> 
{ using type = Opm::SimpleTest<TypeTag>; };

// Set flash solver
template <class TypeTag>
struct FlashSolver<TypeTag, TTag::SimpleTest> {
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using type = Opm::PTFlash<Scalar, FluidSystem>;
};

// Set fluid configuration
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::SimpleTest> 
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::ThreeComponentFluidSystem<Scalar>;
};

// Set the material Law
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::SimpleTest> 
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               //    /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx, TODO
                                               /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                               /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    // define the material law which is parameterized by effective saturations

    using EffMaterialLaw = Opm::NullMaterial<Traits>;
    //using EffMaterialLaw = Opm::BrooksCorey<Traits>;

public:
     using type = EffMaterialLaw;
};

// Write the Newton convergence behavior to disk?
template <class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::SimpleTest> {
static constexpr bool value = false; };

// Enable gravity false
template <class TypeTag>
struct EnableGravity<TypeTag, TTag::SimpleTest> { static constexpr bool value = false;
};

// set the defaults for the problem specific properties
 template <class TypeTag>
 struct Temperature<TypeTag, TTag::SimpleTest> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;//TODO
 };

template <class TypeTag>
struct Initialpressure<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e5;
};

template <class TypeTag>
struct SimulationName<TypeTag, TTag::SimpleTest> {
    static constexpr auto value = "simpletest";
};

// The default for the end time of the simulation
template <class TypeTag>
struct EndTime<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3 * 24. * 60. * 60.;//3 days
};

// convergence control
template <class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    //static constexpr type value = 30;
    static constexpr type value = 1 * 24. * 60. * 60.;
};

template <class TypeTag>
struct LinearSolverTolerance<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

template <class TypeTag>
struct LinearSolverAbsTolerance<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.;
};

template <class TypeTag>
struct NewtonTolerance<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

template <class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60 * 60;
};

template <class TypeTag>
struct NewtonMaxIterations<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 30;
};

template <class TypeTag>
struct NewtonTargetIterations<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 6;
};

// output
template <class TypeTag>
struct VtkWriteFilterVelocities<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWritePotentialGradients<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWriteTotalMassFractions<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWriteTotalMoleFractions<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWriteFugacityCoeffs<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWriteLiquidMoleFractions<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct VtkWriteEquilibriumConstants<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = true;
};

// write restart for every hour
template <class TypeTag>
struct EpisodeLength<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60. * 60.;
};

// mesh grid
template <class TypeTag>
struct Vanguard<TypeTag, TTag::SimpleTest> {
    using type = Opm::StructuredGridVanguard<TypeTag>;
};

//\Note: from the Julia code, the problem is a 1D problem with 3X1 cell.
//\Note: DomainSizeX is 3.0 meters
//\Note: DomainSizeY is 1.0 meters
template <class TypeTag>
struct DomainSizeX<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3; // meter
};

template <class TypeTag>
struct DomainSizeY<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template <class TypeTag>
struct DomainSizeZ<TypeTag, TTag::SimpleTest> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::SimpleTest> { static constexpr int value = 3; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::SimpleTest> { static constexpr int value = 1; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::SimpleTest> { static constexpr int value = 1; };


// compositional, with diffusion
template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::SimpleTest> {
    static constexpr bool value = false;
};

} // namespace Opm::Properties




namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief 3 component simple testproblem with ["CO2", "C1", "C10"]
 *
 */
template <class TypeTag>
class SimpleTest : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;    

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { Comp2Idx = FluidSystem::Comp2Idx };
    enum { Comp1Idx = FluidSystem::Comp1Idx };
    enum { Comp0Idx = FluidSystem::Comp0Idx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + Comp1Idx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableGravity = getPropValue<TypeTag, Properties::EnableGravity>() };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;

public:
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    SimpleTest(Simulator& simulator)
        : ParentType(simulator)
    {
        Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
        simulator.setEpisodeLength(epi_len);
    }

    void initPetrophysics()
    {
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        K_ = this->toDimMatrix_(9.869232667160131e-14);

        porosity_ = 0.1;
    
    }

    template <class Context>
    const DimVector&
    gravity([[maybe_unused]]const Context& context, [[maybe_unused]] unsigned spaceIdx, [[maybe_unused]] unsigned timeIdx) const
    {
        return gravity();
    }

    const DimVector& gravity() const
    {
        return gravity_;
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();
        // initialize fixed parameters; temperature, permeability, porosity
        initPetrophysics();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature, 
                            "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Initialpressure, 
                            "The initial pressure [Pa s] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag,
                             std::string,
                             SimulationName,
                             "The name of the simulation used for the output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EpisodeLength, 
                            "Time interval [s] for episode length");
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName);
        return oss.str();
    }

    // This method must be overridden for the simulator to continue with
    // a new episode. We just start a new episode with the same length as
    // the old one.
    void endEpisode()
    {
        Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
        this->simulator().startNextEpisode(epi_len);
    }

    // only write output when episodes change, aka. report steps, and
    // include the initial timestep too
    bool shouldWriteOutput()
    {
        return this->simulator().episodeWillBeOver() || (this->simulator().timeStepIndex() == -1);
    }

    // we don't care about doing restarts from every fifth timestep, it
    // will just slow us down
    bool shouldWriteRestartFile()
    {
        return false;
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
        Scalar tol = this->model().newtonMethod().tolerance() * 1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageO, storageW;
        this->model().globalPhaseStorage(storageO, oilPhaseIdx);

        // Calculate storage terms
         PrimaryVariables storageL, storageG;
         this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
         this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

         // Write mass balance information for rank 0
        //  if (this->gridView().comm().rank() == 0) {
        //      std::cout << "Storage: liquid=[" << storageL << "]"
        //                << " gas=[" << storageG << "]\n" << std::flush;
        //  }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
        initialFluidState(fs, context, spaceIdx, timeIdx);
        values.assignNaive(fs);
        std::cout << "primary variables for cell " << context.globalSpaceIndex(spaceIdx, timeIdx) << ": " << values << "\n";
    }

    // Constant temperature
    template <class Context>
    Scalar temperature([[maybe_unused]] const Context& context, [[maybe_unused]] unsigned spaceIdx, [[maybe_unused]] unsigned timeIdx) const
    {
        return temperature_;
    }


    // Constant permeability
    template <class Context>
    const DimMatrix& intrinsicPermeability([[maybe_unused]] const Context& context,
                                           [[maybe_unused]] unsigned spaceIdx,
                                           [[maybe_unused]] unsigned timeIdx) const
    {
        return K_;
    }

    // Constant porosity
    template <class Context>
    Scalar porosity([[maybe_unused]] const Context& context, [[maybe_unused]] unsigned spaceIdx, [[maybe_unused]] unsigned timeIdx) const
    {
        int spatialIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        int inj = 0;
        int prod = 2;
    if (spatialIdx == inj || spatialIdx == prod)
            return 1.0;
        else
            return porosity_;
    }

        /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams([[maybe_unused]] const Context& context,
                                               [[maybe_unused]] unsigned spaceIdx,
                                               [[maybe_unused]] unsigned timeIdx) const
    {
        return this->mat_;
    }


    // No flow (introduce fake wells instead)
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& /*context*/,
                  unsigned /*spaceIdx*/,
                  unsigned /*timeIdx*/) const
    { values.setNoFlow(); }

    // No source terms
    template <class Context>
    void source(RateVector& source_rate,
                [[maybe_unused]] const Context& context,
                [[maybe_unused]] unsigned spaceIdx,
                [[maybe_unused]] unsigned timeIdx) const
    {
        source_rate = Scalar(0.0);
    }

private:

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class FluidState, class Context>
    void initialFluidState(FluidState& fs, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
    using Scalar = double;
    using FluidSystem = Opm::ThreeComponentFluidSystem<Scalar>;

    constexpr auto numComponents = FluidSystem::numComponents;
    using Evaluation = Opm::DenseAd::Evaluation<double, numComponents>;
    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;

    // input from Olav
    //z0 = [0.5, 0.3, 0.2]
    //zi = [0.99, 0.01-1e-3, 1e-3]
    //p0 = 75e5
    //T0 = 423.25
    int inj = 0;
    int prod = 2;
    int spatialIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
    ComponentVector comp;
    comp[0] = Evaluation::createVariable(0.5, 1);
    comp[1] = Evaluation::createVariable(0.3, 2);
    comp[2] = 1. - comp[0] - comp[1];
     if (spatialIdx == inj){
         comp[0] = Evaluation::createVariable(0.99, 1);
         comp[1] = Evaluation::createVariable(0.01-1e-3, 2);
         comp[2] = 1. - comp[0] - comp[1]; 
     }
    ComponentVector sat;
    sat[0] = 1.0; sat[1] = 1.0-sat[0];
    // TODO: should we put the derivative against the temperature here?
    const Scalar temp = 423.25;

    // TODO: no capillary pressure for now
    Scalar p0 = 75e5; //CONVERGENCE FAILURE WITH 75

    //\Note, for an AD variable, if we multiply it with 2, the derivative will also be scalced with 2,
    //\Note, so we should not do it.
    if (spatialIdx == inj){
        p0 *= 2.0;
    }
    if (spatialIdx == prod) {
        p0 *= 0.5;
    }
    Evaluation p_init = Evaluation::createVariable(p0, 0);

    fs.setPressure(FluidSystem::oilPhaseIdx, p_init);
    fs.setPressure(FluidSystem::gasPhaseIdx, p_init);

    fs.setMoleFraction(FluidSystem::oilPhaseIdx, FluidSystem::Comp0Idx, comp[0]);
    fs.setMoleFraction(FluidSystem::oilPhaseIdx, FluidSystem::Comp1Idx, comp[1]);
    fs.setMoleFraction(FluidSystem::oilPhaseIdx, FluidSystem::Comp2Idx, comp[2]);

    fs.setMoleFraction(FluidSystem::gasPhaseIdx, FluidSystem::Comp0Idx, comp[0]);
    fs.setMoleFraction(FluidSystem::gasPhaseIdx, FluidSystem::Comp1Idx, comp[1]);
    fs.setMoleFraction(FluidSystem::gasPhaseIdx, FluidSystem::Comp2Idx, comp[2]);

    // It is used here only for calculate the z
    fs.setSaturation(FluidSystem::oilPhaseIdx, sat[0]);
    fs.setSaturation(FluidSystem::gasPhaseIdx, sat[1]);

    fs.setTemperature(temp);

    // ParameterCache paramCache;
    {
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fs, FluidSystem::oilPhaseIdx);
        paramCache.updatePhase(fs, FluidSystem::gasPhaseIdx);
        fs.setDensity(FluidSystem::oilPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::oilPhaseIdx));
        fs.setDensity(FluidSystem::gasPhaseIdx, FluidSystem::density(fs, paramCache, FluidSystem::gasPhaseIdx));
        fs.setViscosity(FluidSystem::oilPhaseIdx, FluidSystem::viscosity(fs, paramCache, FluidSystem::oilPhaseIdx));
        fs.setViscosity(FluidSystem::gasPhaseIdx, FluidSystem::viscosity(fs, paramCache, FluidSystem::gasPhaseIdx));
    }

    // ComponentVector zInit(0.); // TODO; zInit needs to be normalized.
    // {
    //     Scalar sumMoles = 0.0;
    //     for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
    //         for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
    //             Scalar tmp = Opm::getValue(fs.molarity(phaseIdx, compIdx) * fs.saturation(phaseIdx));
    //             zInit[compIdx] += Opm::max(tmp, 1e-8);
    //             sumMoles += tmp;
    //         }
    //     }
    //     zInit /= sumMoles;
    //     // initialize the derivatives
    //     // TODO: the derivative eventually should be from the reservoir flow equations
    //     Evaluation z_last = 1.;
    //     for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
    //         zInit[compIdx] = Evaluation::createVariable(Opm::getValue(zInit[compIdx]), compIdx + 1);
    //         z_last -= zInit[compIdx];
    //     }
    //     zInit[numComponents - 1] = z_last;
    // }

    // TODO: only, p, z need the derivatives.
    const double flash_tolerance = 1.e-12; // just to test the setup in co2-compositional
    //const int flash_verbosity = 1;
    const std::string flash_twophase_method = "newton"; // "ssi"
    //const std::string flash_twophase_method = "ssi";
    // const std::string flash_twophase_method = "ssi+newton";

    // Set initial K and L
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        const Evaluation Ktmp = fs.wilsonK_(compIdx);
        fs.setKvalue(compIdx, Ktmp);
    }

    const Evaluation& Ltmp = -1.0;
    fs.setLvalue(Ltmp);
  
    }

    DimMatrix K_;
    Scalar porosity_;
    Scalar temperature_;
    MaterialLawParams mat_;
    DimVector gravity_;
};
} // namespace Opm

#endif