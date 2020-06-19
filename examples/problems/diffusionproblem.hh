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
 * \copydoc Opm::DiffusionProblem
 */
#ifndef EWOMS_POWER_INJECTION_PROBLEM_HH
#define EWOMS_POWER_INJECTION_PROBLEM_HH

#include <opm/models/ncp/ncpproperties.hh>

#include <opm/models/io/cubegridvanguard.hh>

#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class DiffusionProblem;
}

namespace Opm::Properties {

namespace TTag {

struct DiffusionBaseProblem {};

} // namespace TTag

// Set the grid implementation to be used
template<class TypeTag>
struct Grid<TypeTag, TTag::DiffusionBaseProblem> { using type = Dune::YaspGrid</*dim=*/1>; };

// set the Vanguard property
template<class TypeTag>
struct Vanguard<TypeTag, TTag::DiffusionBaseProblem> { using type = Opm::CubeGridVanguard<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DiffusionBaseProblem> { using type = Opm::DiffusionProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DiffusionBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::H2ON2FluidSystem<Scalar>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::DiffusionBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static_assert(FluidSystem::numPhases == 2,
                  "A fluid system with two phases is required "
                  "for this problem!");

    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>;

public:
    using type = Opm::LinearMaterial<Traits>;
};

// Enable molecular diffusion for this problem
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::DiffusionBaseProblem> { static constexpr bool value = true; };

// Disable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::DiffusionBaseProblem> { static constexpr bool value = false; };

// define the properties specific for the diffusion problem
template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::DiffusionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::DiffusionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::DiffusionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0;
};

template<class TypeTag>
struct CellsX<TypeTag, TTag::DiffusionBaseProblem> { static constexpr int value = 250; };
template<class TypeTag>
struct CellsY<TypeTag, TTag::DiffusionBaseProblem> { static constexpr int value = 1; };
template<class TypeTag>
struct CellsZ<TypeTag, TTag::DiffusionBaseProblem> { static constexpr int value = 1; };

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::DiffusionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e6;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::DiffusionBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1000;
};

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 * \brief 1D problem which is driven by molecular diffusion.
 *
 * The domain is one meter long and completely filled with gas and
 * closed on all boundaries. Its left half exhibits a slightly higher
 * water concentration than the right one. After a while, the
 * concentration of water will be equilibrate due to molecular
 * diffusion.
 */
template <class TypeTag>
class DiffusionProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;

    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    DiffusionProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        FluidSystem::init();

        temperature_ = 273.15 + 20.0;

        materialParams_.finalize();

        K_ = this->toDimMatrix_(1e-12); // [m^2]

        setupInitialFluidStates_();
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return std::string("diffusion_") + Model::name(); }

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
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return 0.35; }

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
     * This problem sets no-flow boundaries everywhere.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context OPM_UNUSED,
                  unsigned spaceIdx OPM_UNUSED,
                  unsigned timeIdx OPM_UNUSED) const
    { values.setNoFlow(); }

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
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (onLeftSide_(pos))
            values.assignNaive(leftInitialFluidState_);
        else
            values.assignNaive(rightInitialFluidState_);
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
    bool onLeftSide_(const GlobalPosition& pos) const
    { return pos[0] < (this->boundingBoxMin()[0] + this->boundingBoxMax()[0]) / 2; }

    void setupInitialFluidStates_()
    {
        // create the initial fluid state for the left half of the domain
        leftInitialFluidState_.setTemperature(temperature_);

        Scalar Sl = 0.0;
        leftInitialFluidState_.setSaturation(liquidPhaseIdx, Sl);
        leftInitialFluidState_.setSaturation(gasPhaseIdx, 1 - Sl);

        Scalar p = 1e5;
        leftInitialFluidState_.setPressure(liquidPhaseIdx, p);
        leftInitialFluidState_.setPressure(gasPhaseIdx, p);

        Scalar xH2O = 0.01;
        leftInitialFluidState_.setMoleFraction(gasPhaseIdx, H2OIdx, xH2O);
        leftInitialFluidState_.setMoleFraction(gasPhaseIdx, N2Idx, 1 - xH2O);

        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        CFRP::solve(leftInitialFluidState_, paramCache, gasPhaseIdx,
                    /*setViscosity=*/false, /*setEnthalpy=*/false);

        // create the initial fluid state for the right half of the domain
        rightInitialFluidState_.assign(leftInitialFluidState_);
        xH2O = 0.0;
        rightInitialFluidState_.setMoleFraction(gasPhaseIdx, H2OIdx, xH2O);
        rightInitialFluidState_.setMoleFraction(gasPhaseIdx, N2Idx, 1 - xH2O);
        CFRP::solve(rightInitialFluidState_, paramCache, gasPhaseIdx,
                    /*setViscosity=*/false, /*setEnthalpy=*/false);
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> leftInitialFluidState_;
    Opm::CompositionalFluidState<Scalar, FluidSystem> rightInitialFluidState_;
    Scalar temperature_;
};

} // namespace Opm

#endif
