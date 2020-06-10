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
 * \copydoc Opm::Tutorial1Problem
 */
#ifndef EWOMS_TUTORIAL1_PROBLEM_HH /*@\label{tutorial1:guardian1}@*/
#define EWOMS_TUTORIAL1_PROBLEM_HH /*@\label{tutorial1:guardian2}@*/

// The numerical model
#include <opm/models/immiscible/immisciblemodel.hh>

// The spatial discretization (VCFV == Vertex-Centered Finite Volumes)
#include <opm/models/discretization/vcfv/vcfvdiscretization.hh>  /*@\label{tutorial1:include-discretization}@*/

// The chemical species that are used
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Lnapl.hpp>

// Headers required for the capillary pressure law
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp> /*@\label{tutorial1:rawLawInclude}@*/
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

// For the DUNE grid
#include <dune/grid/yaspgrid.hh> /*@\label{tutorial1:include-grid-manager}@*/
#include <opm/models/io/cubegridvanguard.hh> /*@\label{tutorial1:include-grid-manager}@*/

// For Dune::FieldMatrix
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Opm {
// forward declaration of the problem class
template <class TypeTag>
class Tutorial1Problem;
}

namespace Opm::Properties {

// Create a new type tag for the problem
// Create new type tags
namespace TTag {
struct Tutorial1Problem { using InheritsFrom = std::tuple<ImmiscibleTwoPhaseModel>; };
} // end namespace TTag

// Select the vertex centered finite volume method as spatial discretization
template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Tutorial1Problem>
{ using type = TTag::VcfvDiscretization; }; /*@\label{tutorial1:set-spatial-discretization}@*/

// Set the "Problem" property
template<class TypeTag>
struct Problem<TypeTag, TTag::Tutorial1Problem>
{ using type = Opm::Tutorial1Problem<TypeTag>; }; /*@\label{tutorial1:set-problem}@*/

// Set grid and the grid manager to be used
template<class TypeTag>
struct Grid<TypeTag, TTag::Tutorial1Problem> { using type = Dune::YaspGrid</*dim=*/2>; }; /*@\label{tutorial1:set-grid}@*/
template<class TypeTag>
struct Vanguard<TypeTag, TTag::Tutorial1Problem> { using type = Opm::CubeGridVanguard<TypeTag>; }; /*@\label{tutorial1:set-grid-manager}@*/

// Set the wetting phase /*@\label{tutorial1:2p-system-start}@*/
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::Tutorial1Problem> /*@\label{tutorial1:wettingPhase}@*/
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::Tutorial1Problem> /*@\label{tutorial1:nonwettingPhase}@*/
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Opm::LiquidPhase<Scalar, Opm::LNAPL<Scalar> >;
}; /*@\label{tutorial1:2p-system-end}@*/

// Set the material law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Tutorial1Problem>
{
private:
    // create a class holding the necessary information for a
    // two-phase capillary pressure law
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar, wettingPhaseIdx, nonWettingPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using RawMaterialLaw = Opm::RegularizedBrooksCorey<Traits>; /*@\label{tutorial1:rawlaw}@*/

public:
    // Convert absolute saturations into effective ones before passing
    // it to the base capillary pressure law
    using type = Opm::EffToAbsLaw<RawMaterialLaw>; /*@\label{tutorial1:eff2abs}@*/
};

// Disable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::Tutorial1Problem> { static constexpr bool value = false; }; /*@\label{tutorial1:gravity}@*/

// define how long the simulation should run [s]
template<class TypeTag>
struct EndTime<TypeTag, TTag::Tutorial1Problem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100e3;
}; /*@\label{tutorial1:default-params-begin}@*/

// define the size of the initial time step [s]
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Tutorial1Problem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 125.0;
};

// define the physical size of the problem's domain [m]
template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::Tutorial1Problem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 300.0;
}; /*@\label{tutorial1:grid-default-params-begin}@*/
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::Tutorial1Problem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60.0;
};
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::Tutorial1Problem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

// // define the number of cells used for discretizing the physical domain
template<class TypeTag>
struct CellsX<TypeTag, TTag::Tutorial1Problem> { static constexpr int value = 100; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::Tutorial1Problem> { static constexpr int value = 1; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::Tutorial1Problem> { static constexpr int value = 1; }; /*@\label{tutorial1:default-params-end}@*/

} // namespace Opm::Properties

namespace Opm {
//! Tutorial problem using the "immiscible" model.
template <class TypeTag>
class Tutorial1Problem
    : public GetPropType<TypeTag, Properties::BaseProblem> /*@\label{tutorial1:def-problem}@*/
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    // Grid dimension
    enum { dimWorld = GridView::dimensionworld };

    // The type of the intrinsic permeability tensor
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    // eWoms specific types are specified via the property system
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>; /*@\label{tutorial1:matLawObjectType}@*/

    // phase indices
    enum { numPhases = FluidSystem::numPhases };
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    // Indices of the conservation equations
    enum { contiWettingEqIdx = Indices::conti0EqIdx + wettingPhaseIdx };
    enum { contiNonWettingEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx };

public:
    //! The constructor of the problem. This only _allocates_ the memory required by the
    //! problem. The constructor is supposed to _never ever_ throw an exception.
    Tutorial1Problem(Simulator& simulator)
        : ParentType(simulator)
        , eps_(3e-6)
    { }

    //! This method initializes the data structures allocated by the problem
    //! constructor. In contrast to the constructor, exceptions thrown from within this
    //! method won't lead to segmentation faults.
    void finishInit()
    {
        ParentType::finishInit();

        // Use an isotropic and homogeneous intrinsic permeability
        K_ = this->toDimMatrix_(1e-7);

        // Parameters of the Brooks-Corey law
        materialParams_.setEntryPressure(500.0 /*Pa*/); /*@\label{tutorial1:setLawParams}@*/
        materialParams_.setLambda(2); // shape parameter

        // Set the residual saturations
        materialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        materialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

        // wrap up the initialization of the material law's parameters
        materialParams_.finalize();
    }

    //! Specifies the problem name. This is used for files generated by the simulation.
    std::string name() const
    { return "tutorial1"; }

    //! Returns the temperature at a given position.
    template <class Context>
    Scalar temperature(const Context& /*context*/,
                       unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return 283.15; }

    //! Returns the intrinsic permeability tensor [m^2] at a position.
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& /*context*/, /*@\label{tutorial1:permeability}@*/
                                           unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return K_; }

    //! Defines the porosity [-] of the medium at a given position
    template <class Context>
    Scalar porosity(const Context& /*context*/,
                    unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const /*@\label{tutorial1:porosity}@*/
    { return 0.2; }

    //! Returns the parameter object for the material law at a given position
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& /*context*/, /*@\label{tutorial1:matLawParams}@*/
                                               unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return materialParams_; }

    //! Evaluates the boundary conditions.
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (pos[0] < eps_) {
            // Free-flow conditions on left boundary
            const auto& materialParams = this->materialLawParams(context, spaceIdx, timeIdx);

            Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            Scalar Sw = 1.0;
            fs.setSaturation(wettingPhaseIdx, Sw);
            fs.setSaturation(nonWettingPhaseIdx, 1.0 - Sw);
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, materialParams, fs);
            fs.setPressure(wettingPhaseIdx, 200e3);
            fs.setPressure(nonWettingPhaseIdx, 200e3 + pC[nonWettingPhaseIdx] - pC[nonWettingPhaseIdx]);

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
                fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
            }

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (pos[0] > this->boundingBoxMax()[0] - eps_) {
            // forced outflow at the right boundary
            RateVector massRate(0.0);

            massRate[contiWettingEqIdx] = 0.0;  // [kg / (s m^2)]
            massRate[contiNonWettingEqIdx] = 3e-2; // [kg / (s m^2)]

            values.setMassRate(massRate);
        }
        else // no flow at the remaining boundaries
            values.setNoFlow();
    }

    //! Evaluates the source term for all conserved quantities at a given
    //! position of the domain [kg/(m^3 * s)]. Positive values mean that
    //! mass is created.
    template <class Context>
    void source(RateVector& sourceRate, const Context& /*context*/,
                unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    {
        sourceRate[contiWettingEqIdx] = 0.0;
        sourceRate[contiNonWettingEqIdx] = 0.0;
    }

    //! Evaluates the initial value at a given position in the domain.
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context,
                 unsigned spaceIdx, unsigned timeIdx) const
    {
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        // the domain is initially fully saturated by LNAPL
        Scalar Sw = 0.0;
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1.0 - Sw);

        // the temperature is given by the temperature() method
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        // set pressure of the wetting phase to 200 kPa = 2 bar
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, materialLawParams(context, spaceIdx, timeIdx),
                                        fs);
        fs.setPressure(wettingPhaseIdx, 200e3);
        fs.setPressure(nonWettingPhaseIdx, 200e3 + pC[nonWettingPhaseIdx] - pC[nonWettingPhaseIdx]);

        values.assignNaive(fs);
    }

private:
    DimMatrix K_;
    // Object that holds the parameters of required by the capillary pressure law.
    MaterialLawParams materialParams_; /*@\label{tutorial1:matParamsObject}@*/

    // small epsilon value
    Scalar eps_;
};
} // namespace Opm

#endif

