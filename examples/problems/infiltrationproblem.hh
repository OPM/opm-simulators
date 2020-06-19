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
 * \copydoc Opm::InfiltrationProblem
 */
#ifndef EWOMS_INFILTRATION_PROBLEM_HH
#define EWOMS_INFILTRATION_PROBLEM_HH

#include <opm/models/pvs/pvsproperties.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/ThreePhaseParkerVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class InfiltrationProblem;
}

namespace Opm::Properties {

namespace TTag {
struct InfiltrationBaseProblem {};
}

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationBaseProblem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationBaseProblem> { using type = Opm::InfiltrationProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationBaseProblem>
{ using type = Opm::H2OAirMesityleneFluidSystem<GetPropType<TypeTag, Properties::Scalar>>; };

// Enable gravity?
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::InfiltrationBaseProblem> { static constexpr bool value = true; };

// Write newton convergence?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::InfiltrationBaseProblem> { static constexpr bool value = false; };

// -1 backward differences, 0: central differences, +1: forward differences
template<class TypeTag>
struct NumericDifferenceMethod<TypeTag, TTag::InfiltrationBaseProblem> { static constexpr int value = 1; };

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::InfiltrationBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using Traits= Opm::ThreePhaseMaterialTraits<
        Scalar,
        /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
        /*nonWettingPhaseIdx=*/FluidSystem::naplPhaseIdx,
        /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx>;

public:
    using type = Opm::ThreePhaseParkerVanGenuchten<Traits>;
};

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::InfiltrationBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 6e3;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::InfiltrationBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60;
};

// The default DGF file to load
template<class TypeTag>
struct GridFile<TypeTag, TTag::InfiltrationBaseProblem>
{ static constexpr auto value = "./data/infiltration_50x3.dgf"; };

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 * \brief Isothermal NAPL infiltration problem where LNAPL
 *        contaminates the unsaturated and the saturated groundwater
 *        zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep,
 * where the lower part represents a slightly inclined groundwater
 * table, and the upper part is the vadose zone.  A LNAPL (Non-Aqueous
 * Phase Liquid which is lighter than water) infiltrates (modelled
 * with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and
 * migrates on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually
 * trapped immobile NAPL, which can in the following dissolve and
 * evaporate slowly, and eventually be transported by advection and
 * diffusion.
 *
 * Left and right boundaries are constant hydraulic head boundaries
 * (Dirichlet), Top and bottom are Neumann boundaries, all no-flow
 * except for the small infiltration zone in the upper left part.
 */
template <class TypeTag>
class InfiltrationProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    // copy some indices for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum {
        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,

        // number of phases/components
        numPhases = FluidSystem::numPhases,

        // component indices
        NAPLIdx = FluidSystem::NAPLIdx,
        H2OIdx = FluidSystem::H2OIdx,
        airIdx = FluidSystem::airIdx,

        // phase indices
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        naplPhaseIdx = FluidSystem::naplPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    InfiltrationProblem(Simulator& simulator)
        : ParentType(simulator)
        , eps_(1e-6)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(/*tempMin=*/temperature_ - 1,
                          /*tempMax=*/temperature_ + 1,
                          /*nTemp=*/3,
                          /*pressMin=*/0.8 * 1e5,
                          /*pressMax=*/3 * 1e5,
                          /*nPress=*/200);

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-11);
        coarseK_ = this->toDimMatrix_(1e-11);

        // porosities
        porosity_ = 0.40;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSwrx(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the three-phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgN(4.);
        materialParams_.setkrRegardsSnr(false);

        materialParams_.finalize();
        materialParams_.checkDefined();
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::shouldWriteRestartFile
     *
     * This problem writes a restart file after every time step.
     */
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "infiltration_" << Model::name();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

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

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix&
    intrinsicPermeability(const Context& context,
                          unsigned spaceIdx,
                          unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return porosity_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& context OPM_UNUSED,
                      unsigned spaceIdx OPM_UNUSED,
                      unsigned timeIdx OPM_UNUSED) const
    { return materialParams_; }

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

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector molarRate(0.0);
            molarRate[conti0EqIdx + NAPLIdx] = -0.001;

            values.setMolarRate(molarRate);
            Opm::Valgrind::CheckDefined(values);
        }
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
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        Opm::Valgrind::CheckDefined(values);
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

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onUpperBoundary_(pos) && 50 < pos[0] && pos[0] < 75; }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState& fs, const Context& context,
                            unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition pos = context.pos(spaceIdx, timeIdx);
        Scalar y = pos[1];
        Scalar x = pos[0];

        Scalar densityW = 1000.0;
        Scalar pc = 9.81 * densityW * (y - (5 - 5e-4 * x));
        if (pc < 0.0)
            pc = 0.0;

        // set pressures
        const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
        Scalar Sw = matParams.Swr();
        Scalar Swr = matParams.Swr();
        Scalar Sgr = matParams.Sgr();
        if (Sw < Swr)
            Sw = Swr;
        if (Sw > 1 - Sgr)
            Sw = 1 - Sgr;
        Scalar Sg = 1 - Sw;

        Opm::Valgrind::CheckDefined(Sw);
        Opm::Valgrind::CheckDefined(Sg);

        fs.setSaturation(waterPhaseIdx, Sw);
        fs.setSaturation(gasPhaseIdx, Sg);
        fs.setSaturation(naplPhaseIdx, 0);

        // set temperature of all phases
        fs.setTemperature(temperature_);

        // compute pressures
        Scalar pcAll[numPhases];
        Scalar pg = 1e5;
        if (onLeftBoundary_(pos))
            pg += 10e3;
        MaterialLaw::capillaryPressures(pcAll, matParams, fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fs.setPressure(phaseIdx, pg + (pcAll[phaseIdx] - pcAll[gasPhaseIdx]));

        // set composition of gas phase
        fs.setMoleFraction(gasPhaseIdx, H2OIdx, 1e-6);
        fs.setMoleFraction(gasPhaseIdx, airIdx,
                           1 - fs.moleFraction(gasPhaseIdx, H2OIdx));
        fs.setMoleFraction(gasPhaseIdx, NAPLIdx, 0);

        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        CFRP::solve(fs, paramCache, gasPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/false);

        fs.setMoleFraction(waterPhaseIdx, H2OIdx,
                           1 - fs.moleFraction(waterPhaseIdx, H2OIdx));
    }

    bool isFineMaterial_(const GlobalPosition& pos) const
    {  return 70. <= pos[0] && pos[0] <= 85. && 7.0 <= pos[1] && pos[1] <= 7.50; }

    DimMatrix fineK_;
    DimMatrix coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Opm

#endif
