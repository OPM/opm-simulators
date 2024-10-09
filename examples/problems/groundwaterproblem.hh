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
 * \copydoc Opm::GroundWaterProblem
 */
#ifndef EWOMS_GROUND_WATER_PROBLEM_HH
#define EWOMS_GROUND_WATER_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>

#include <opm/models/common/multiphasebaseparameters.hh>

#include <opm/models/immiscible/immiscibleproperties.hh>

#include <opm/simulators/linalg/parallelistlbackend.hh>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class GroundWaterProblem;
}

namespace Opm::Properties {

namespace TTag {
struct GroundWaterBaseProblem {};
}

template<class TypeTag>
struct Fluid<TypeTag, TTag::GroundWaterBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::GroundWaterBaseProblem> { using type = Dune::YaspGrid<2>; };
// struct Grid<TypeTag, TTag::GroundWaterBaseProblem> { using type = Dune::SGrid<2, 2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::GroundWaterBaseProblem>
{ using type = Opm::GroundWaterProblem<TypeTag>; };

// Use the conjugated gradient linear solver with the default preconditioner (i.e.,
// ILU-0) from dune-istl
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::GroundWaterBaseProblem> { using type = TTag::ParallelIstlLinearSolver; };

template<class TypeTag>
struct LinearSolverWrapper<TypeTag, TTag::GroundWaterBaseProblem>
{ using type = Opm::Linear::SolverWrapperConjugatedGradients<TypeTag>; };

} // namespace Opm::Properties

namespace Opm::Parameters {

template<class Scalar>
struct LensLowerLeftX { static constexpr Scalar value = 0.25; };

template<class Scalar>
struct LensLowerLeftY { static constexpr Scalar value = 0.25; };

template<class Scalar>
struct LensLowerLeftZ { static constexpr Scalar value = 0.25; };

template<class Scalar>
struct LensUpperRightX { static constexpr Scalar value = 0.75; };

template<class Scalar>
struct LensUpperRightY { static constexpr Scalar value = 0.75; };

template<class Scalar>
struct LensUpperRightZ { static constexpr Scalar value = 0.75; };

template<class Scalar>
struct Permeability { static constexpr Scalar value = 1e-10; };

template<class Scalar>
struct PermeabilityLens { static constexpr Scalar value = 1e-12; };

} // namespace Opm::Parameters

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Test for the immisicible VCVF discretization with only a single phase
 *
 * This problem is inspired by groundwater flow. Don't expect it to be
 * realistic, though: For two dimensions, the domain size is 1m times
 * 1m. On the left and right of the domain, no-flow boundaries are
 * used, while at the top and bottom free flow boundaries with a
 * pressure of 2 bar and 1 bar are used. The center of the domain is
 * occupied by a rectangular lens of lower permeability.
 */
template <class TypeTag>
class GroundWaterProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // copy some indices for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum {
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressure0Idx = Indices::pressure0Idx
    };

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    GroundWaterProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1.0e-3;

        lensLowerLeft_[0] = Parameters::Get<Parameters::LensLowerLeftX<Scalar>>();
        if (dim > 1)
            lensLowerLeft_[1] = Parameters::Get<Parameters::LensLowerLeftY<Scalar>>();
        if (dim > 2)
            lensLowerLeft_[2] = Parameters::Get<Parameters::LensLowerLeftY<Scalar>>();

        lensUpperRight_[0] = Parameters::Get<Parameters::LensUpperRightX<Scalar>>();
        if (dim > 1)
            lensUpperRight_[1] = Parameters::Get<Parameters::LensUpperRightY<Scalar>>();
        if (dim > 2)
            lensUpperRight_[2] = Parameters::Get<Parameters::LensUpperRightY<Scalar>>();

        intrinsicPerm_ = this->toDimMatrix_(Parameters::Get<Parameters::Permeability<Scalar>>());
        intrinsicPermLens_ = this->toDimMatrix_(Parameters::Get<Parameters::PermeabilityLens<Scalar>>());
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        Parameters::Register<Parameters::LensLowerLeftX<Scalar>>
            ("The x-coordinate of the lens' lower-left corner [m].");
        Parameters::Register<Parameters::LensUpperRightX<Scalar>>
            ("The x-coordinate of the lens' upper-right corner [m].");

        if (dimWorld > 1) {
            Parameters::Register<Parameters::LensLowerLeftY<Scalar>>
                ("The y-coordinate of the lens' lower-left corner [m].");
            Parameters::Register<Parameters::LensUpperRightY<Scalar>>
                ("The y-coordinate of the lens' upper-right corner [m].");
        }

        if (dimWorld > 2) {
            Parameters::Register<Parameters::LensLowerLeftZ<Scalar>>
                ("The z-coordinate of the lens' lower-left corner [m].");
            Parameters::Register<Parameters::LensUpperRightZ<Scalar>>
                ("The z-coordinate of the lens' upper-right corner [m].");
        }

        Parameters::Register<Parameters::Permeability<Scalar>>
            ("The intrinsic permeability [m^2] of the ambient material.");
        Parameters::Register<Parameters::PermeabilityLens<Scalar>>
            ("The intrinsic permeability [m^2] of the lens.");

        Parameters::SetDefault<Parameters::GridFile>("./data/groundwater_2d.dgf");
        Parameters::SetDefault<Parameters::EndTime<Scalar>>(1.0);
        Parameters::SetDefault<Parameters::InitialTimeStepSize<Scalar>>(1.0);
        Parameters::SetDefault<Parameters::EnableGravity>(true);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "groundwater_" << Model::name();
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
    Scalar temperature(const Context& /*context*/,
                       unsigned /*spaceIdx*/,
                       unsigned /*timeIdx*/) const
    { return 273.15 + 10; } // 10C

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& /*context*/,
                    unsigned /*spaceIdx*/,
                    unsigned /*timeIdx*/) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        if (isInLens_(context.pos(spaceIdx, timeIdx)))
            return intrinsicPermLens_;
        else
            return intrinsicPerm_;
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
        const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);

        if (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos)) {
            Scalar pressure;
            Scalar T = temperature(context, spaceIdx, timeIdx);
            if (onLowerBoundary_(globalPos))
                pressure = 2e5;
            else // on upper boundary
                pressure = 1e5;

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                                      /*storeEnthalpy=*/false> fs;
            fs.setSaturation(/*phaseIdx=*/0, 1.0);
            fs.setPressure(/*phaseIdx=*/0, pressure);
            fs.setTemperature(T);

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
                fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
            }

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else {
            // no flow boundary
            values.setNoFlow();
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
    void initial(PrimaryVariables& values,
                 const Context& /*context*/,
                 unsigned /*spaceIdx*/,
                 unsigned /*timeIdx*/) const
    {
        // const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);
        values[pressure0Idx] = 1.0e+5; // + 9.81*1.23*(20-globalPos[dim-1]);
    }

    /*!
     * \copydoc FvBaseProblem::source
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& /*context*/,
                unsigned /*spaceIdx*/,
                unsigned /*timeIdx*/) const
    { rate = Scalar(0.0); }

    //! \}

private:
    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[dim - 1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[dim - 1] > this->boundingBoxMax()[dim - 1] - eps_; }

    bool isInLens_(const GlobalPosition& pos) const
    {
        return lensLowerLeft_[0] <= pos[0] && pos[0] <= lensUpperRight_[0]
               && lensLowerLeft_[1] <= pos[1] && pos[1] <= lensUpperRight_[1];
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    DimMatrix intrinsicPerm_;
    DimMatrix intrinsicPermLens_;

    Scalar eps_;
};
} // namespace Opm

#endif
