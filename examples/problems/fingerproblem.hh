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
 * \copydoc Opm::FingerProblem
 */
#ifndef EWOMS_FINGER_PROBLEM_HH
#define EWOMS_FINGER_PROBLEM_HH

#include <opm/models/io/structuredgridvanguard.hh>

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/ParkerLenhard.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Air.hpp>

#include <opm/models/immiscible/immiscibleproperties.hh>
#include <opm/models/discretization/common/restrictprolong.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <vector>
#include <string>

namespace Opm {
template <class TypeTag>
class FingerProblem;

} // namespace Opm

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct FingerBaseProblem { using InheritsFrom = std::tuple<StructuredGridVanguard>; };
} // end namespace TTag

#if HAVE_DUNE_ALUGRID
// use dune-alugrid if available
template<class TypeTag>
struct Grid<TypeTag, TTag::FingerBaseProblem>
{ using type = Dune::ALUGrid</*dim=*/2,
                             /*dimWorld=*/2,
                             Dune::cube,
                             Dune::nonconforming>; };
#endif

// declare the properties used by the finger problem
template<class TypeTag, class MyTypeTag>
struct InitialWaterSaturation { using type = UndefinedProperty; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FingerBaseProblem> { using type = Opm::FingerProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::FingerBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::FingerBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::GasPhase<Scalar, Opm::Air<Scalar> >;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FingerBaseProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>;

    // use the parker-lenhard hysteresis law
    using ParkerLenhard = Opm::ParkerLenhard<Traits>;
    using type = ParkerLenhard;
};

// Write the solutions of individual newton iterations?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::FingerBaseProblem> { static constexpr bool value = false; };

// Use forward differences instead of central differences
template<class TypeTag>
struct NumericDifferenceMethod<TypeTag, TTag::FingerBaseProblem> { static constexpr int value = +1; };

// Enable constraints
template<class TypeTag>
struct EnableConstraints<TypeTag, TTag::FingerBaseProblem> { static constexpr int value = true; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::FingerBaseProblem> { static constexpr bool value = true; };

// define the properties specific for the finger problem
template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1;
};
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.3;
};
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.1;
};

template<class TypeTag>
struct InitialWaterSaturation<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.01;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::FingerBaseProblem> { static constexpr int value = 20; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::FingerBaseProblem> { static constexpr int value = 70; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::FingerBaseProblem> { static constexpr int value = 1; };

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 215;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::FingerBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Two-phase problem featuring some gravity-driven saturation
 *        fingers.
 *
 * The domain of this problem is sized 10cm times 1m and is initially
 * dry. Water is then injected at three locations on the top of the
 * domain which leads to gravity fingering. The boundary conditions
 * used are no-flow for the left and right and top of the domain and
 * free-flow at the bottom. This problem uses the Parker-Lenhard
 * hystersis model which might lead to non-monotonic saturation in the
 * fingers if the right material parameters is chosen and the spatial
 * discretization is fine enough.
 */
template <class TypeTag>
class FingerProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    //!\cond SKIP_THIS
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiWettingEqIdx = Indices::conti0EqIdx + wettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil> ;
    enum { codim = Stencil::Entity::codimension };
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;

    using ParkerLenhard = typename GetProp<TypeTag, Properties::MaterialLaw>::ParkerLenhard;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Grid = typename GridView :: Grid;

    using MaterialLawParamsContainer = Dune::PersistentContainer< Grid, std::shared_ptr< MaterialLawParams > >  ;
    //!\endcond

public:
    using RestrictProlongOperator = CopyRestrictProlong< Grid, MaterialLawParamsContainer >;

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    FingerProblem(Simulator& simulator)
        : ParentType(simulator),
          materialParams_( simulator.vanguard().grid(), codim )
    {
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \brief \copydoc FvBaseProblem::restrictProlongOperator
     */
    RestrictProlongOperator restrictProlongOperator()
    {
        return RestrictProlongOperator( materialParams_ );
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return
            std::string("finger") +
            "_" + Model::name() +
            "_" + Model::discretizationName() +
            (this->model().enableGridAdaptation()?"_adaptive":"");
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialWaterSaturation,
                             "The initial saturation in the domain [] of the "
                             "wetting phase");
    }

    /*!
     * \copydoc FvBaseProblem::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;

        temperature_ = 273.15 + 20; // -> 20°C

        FluidSystem::init();

        // parameters for the Van Genuchten law of the main imbibition
        // and the main drainage curves.
        micParams_.setVgAlpha(0.0037);
        micParams_.setVgN(4.7);
        micParams_.finalize();

        mdcParams_.setVgAlpha(0.0037);
        mdcParams_.setVgN(4.7);
        mdcParams_.finalize();

        // initialize the material parameter objects of the individual
        // finite volumes, resize will resize the container to the number of elements
        materialParams_.resize();

        for (auto it = materialParams_.begin(),
                 end = materialParams_.end(); it != end; ++it ) {
            std::shared_ptr< MaterialLawParams >& materialParams = *it ;
            if( ! materialParams )
            {
                materialParams.reset( new MaterialLawParams() );
                materialParams->setMicParams(&micParams_);
                materialParams->setMdcParams(&mdcParams_);
                materialParams->setSwr(0.0);
                materialParams->setSnr(0.1);
                materialParams->finalize();
                ParkerLenhard::reset(*materialParams);
            }
        }

        K_ = this->toDimMatrix_(4.6e-10);

        setupInitialFluidState_();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
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

        // update the history of the hysteresis law
        ElementContext elemCtx(this->simulator());

        auto elemIt = this->gridView().template begin<0>();
        const auto& elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elemCtx.updateAll( elem );
            size_t numDofs = elemCtx.numDof(/*timeIdx=*/0);
            for (unsigned scvIdx = 0; scvIdx < numDofs; ++scvIdx)
            {
                MaterialLawParams& materialParam = materialLawParams( elemCtx, scvIdx, /*timeIdx=*/0 );
                const auto& fs = elemCtx.intensiveQuantities(scvIdx, /*timeIdx=*/0).fluidState();
                ParkerLenhard::update(materialParam, fs);
            }
        }
    }

    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return K_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    MaterialLawParams& materialLawParams(const Context& context,
                                         unsigned spaceIdx, unsigned timeIdx)
    {
        const auto& entity = context.stencil(timeIdx).entity(spaceIdx);
        assert(materialParams_[entity]);
        return *materialParams_[entity];
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& entity = context.stencil(timeIdx).entity( spaceIdx );
        assert(materialParams_[entity]);
        return *materialParams_[entity];
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
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos) || onLowerBoundary_(pos))
            values.setNoFlow();
        else {
            assert(onUpperBoundary_(pos));

            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
        }

        // override the value for the liquid phase by forced
        // imbibition of water on inlet boundary segments
        if (onInlet_(pos)) {
            values[contiWettingEqIdx] = -0.001; // [kg/(m^2 s)]
        }
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
    void initial(PrimaryVariables& values, const Context& /*context*/, unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
    }

    /*!
     * \copydoc FvBaseProblem::constraints
     */
    template <class Context>
    void constraints(Constraints& constraints, const Context& context,
                     unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos) && !onInlet_(pos)) {
            constraints.setActive(true);
            constraints.assignNaive(initialFluidState_);
        }
        else if (onLowerBoundary_(pos)) {
            constraints.setActive(true);
            constraints.assignNaive(initialFluidState_);
        }
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate, const Context& /*context*/,
                unsigned /*spaceIdx*/, unsigned /*timeIdx*/) const
    { rate = Scalar(0.0); }
    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    {
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar lambda = (this->boundingBoxMax()[0] - pos[0]) / width;

        if (!onUpperBoundary_(pos))
            return false;

        Scalar xInject[] = { 0.25, 0.75 };
        Scalar injectLen[] = { 0.1, 0.1 };
        for (unsigned i = 0; i < sizeof(xInject) / sizeof(Scalar); ++i) {
            if (xInject[i] - injectLen[i] / 2 < lambda
                && lambda < xInject[i] + injectLen[i] / 2)
                return true;
        }
        return false;
    }

    void setupInitialFluidState_()
    {
        auto& fs = initialFluidState_;
        fs.setPressure(wettingPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = EWOMS_GET_PARAM(TypeTag, Scalar, InitialWaterSaturation);
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        // set the absolute pressures
        Scalar pn = 1e5;
        fs.setPressure(nonWettingPhaseIdx, pn);
        fs.setPressure(wettingPhaseIdx, pn);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
            fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
        }

    }

    DimMatrix K_;

    typename MaterialLawParams::VanGenuchtenParams micParams_;
    typename MaterialLawParams::VanGenuchtenParams mdcParams_;

    MaterialLawParamsContainer materialParams_;

    Opm::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;

    Scalar temperature_;
    Scalar eps_;
};

} // namespace Opm

#endif
