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
 * \copydoc Opm::PowerInjectionProblem
 */
#ifndef EWOMS_POWER_INJECTION_PROBLEM_HH
#define EWOMS_POWER_INJECTION_PROBLEM_HH

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/io/cubegridvanguard.hh>

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <type_traits>
#include <iostream>

namespace Opm {
template <class TypeTag>
class PowerInjectionProblem;
}

namespace Opm::Properties {

namespace TTag {
struct PowerInjectionBaseProblem {};
}

// Set the grid implementation to be used
template<class TypeTag>
struct Grid<TypeTag, TTag::PowerInjectionBaseProblem> { using type = Dune::YaspGrid</*dim=*/1>; };

// set the Vanguard property
template<class TypeTag>
struct Vanguard<TypeTag, TTag::PowerInjectionBaseProblem> { using type = Opm::CubeGridVanguard<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PowerInjectionBaseProblem> { using type = Opm::PowerInjectionProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::PowerInjectionBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::PowerInjectionBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::GasPhase<Scalar, Opm::Air<Scalar> >;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::PowerInjectionBaseProblem>
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
    using EffectiveLaw = Opm::RegularizedVanGenuchten<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::EffToAbsLaw<EffectiveLaw>;
};

// Write out the filter velocities for this problem
template<class TypeTag>
struct VtkWriteFilterVelocities<TypeTag, TTag::PowerInjectionBaseProblem> { static constexpr bool value = true; };

// Disable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::PowerInjectionBaseProblem> { static constexpr bool value = false; };

// define the properties specific for the power injection problem
template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::PowerInjectionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100.0;
};
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::PowerInjectionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::PowerInjectionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::PowerInjectionBaseProblem> { static constexpr int value = 250; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::PowerInjectionBaseProblem> { static constexpr int value = 1; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::PowerInjectionBaseProblem> { static constexpr int value = 1; };

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::PowerInjectionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 100;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::PowerInjectionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 * \brief 1D Problem with very fast injection of gas on the left.
 *
 * The velocity model is chosen in the .cc file in this problem. The
 * spatial parameters are inspired by the ones given by
 *
 * V. Jambhekar: "Forchheimer Porous-media Flow models -- Numerical
 * Investigation and Comparison with Experimental Data", Master's
 * Thesis at Institute for Modelling Hydraulic and Environmental
 * Systems, University of Stuttgart, 2011
 */
template <class TypeTag>
class PowerInjectionProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    PowerInjectionProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 26.6;

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(0.00045);
        materialParams_.setVgN(7.3);
        materialParams_.finalize();

        K_ = this->toDimMatrix_(5.73e-08); // [m^2]

        setupInitialFluidState_();
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
        oss << "powerinjection_";
        if (std::is_same<GetPropType<TypeTag, Properties::FluxModule>,
                         Opm::DarcyFluxModule<TypeTag> >::value)
            oss << "darcy";
        else
            oss << "forchheimer";

        if (std::is_same<GetPropType<TypeTag, Properties::LocalLinearizerSplice>,
                         Properties::TTag::AutoDiffLocalLinearizer>::value)
            oss << "_" << "ad";
        else
            oss << "_" << "fd";

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
    //! \}

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
    { return K_; }

    /*!
     * \copydoc ForchheimerBaseProblem::ergunCoefficient
     */
    template <class Context>
    Scalar ergunCoefficient(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    { return 0.3866; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return 0.558; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams&
    materialLawParams(const Context& context OPM_UNUSED,
                      unsigned spaceIdx OPM_UNUSED,
                      unsigned timeIdx OPM_UNUSED) const
    { return materialParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * This problem sets a very high injection rate of nitrogen on the
     * left and a free-flow boundary on the right.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiNEqIdx] = -1.00; // kg / (m^2 * s)

            // impose a forced flow boundary
            values.setMassRate(massRate);
        }
        else if (onRightBoundary_(pos))
            // free flow boundary with initial condition on the right
            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
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
                 const Context& context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
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
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    void setupInitialFluidState_()
    {
        initialFluidState_.setTemperature(temperature_);

        Scalar Sw = 1.0;
        initialFluidState_.setSaturation(wettingPhaseIdx, Sw);
        initialFluidState_.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        Scalar p = 1e5;
        initialFluidState_.setPressure(wettingPhaseIdx, p);
        initialFluidState_.setPressure(nonWettingPhaseIdx, p);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(initialFluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            initialFluidState_.setDensity(phaseIdx,
                                          FluidSystem::density(initialFluidState_, paramCache, phaseIdx));
            initialFluidState_.setViscosity(phaseIdx,
                                            FluidSystem::viscosity(initialFluidState_, paramCache, phaseIdx));
        }
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Opm::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;
    Scalar temperature_;
    Scalar eps_;
};

} // namespace Opm

#endif
