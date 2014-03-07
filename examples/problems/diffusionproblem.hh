/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::DiffusionProblem
 */
#ifndef EWOMS_POWER_INJECTION_PROBLEM_HH
#define EWOMS_POWER_INJECTION_PROBLEM_HH

#include <ewoms/models/ncp/ncpproperties.hh>
#include <ewoms/io/cubegridcreator.hh>

#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class DiffusionProblem;
}

namespace Opm {
namespace Properties {

NEW_TYPE_TAG(DiffusionBaseProblem);

// Set the grid implementation to be used
SET_TYPE_PROP(DiffusionBaseProblem, Grid, Dune::YaspGrid</*dim=*/1>);

// set the GridCreator property
SET_TYPE_PROP(DiffusionBaseProblem, GridCreator, Ewoms::CubeGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(DiffusionBaseProblem, Problem, Ewoms::DiffusionProblem<TypeTag>);

// Set the fluid system
SET_PROP(DiffusionBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::FluidSystems::H2ON2<Scalar> type;
};

// Set the material Law
SET_PROP(DiffusionBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static_assert(FluidSystem::numPhases == 2,
                  "A fluid system with two phases is required "
                  "for this problem!");

    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::lPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gPhaseIdx>
    Traits;

public:
    typedef Opm::LinearMaterial<Traits> type;
};

// Enable molecular diffusion for this problem
SET_BOOL_PROP(DiffusionBaseProblem, EnableDiffusion, true);

// Disable gravity
SET_BOOL_PROP(DiffusionBaseProblem, EnableGravity, false);

// define the properties specific for the diffusion problem
SET_SCALAR_PROP(DiffusionBaseProblem, DomainSizeX, 1.0);
SET_SCALAR_PROP(DiffusionBaseProblem, DomainSizeY, 1.0);
SET_SCALAR_PROP(DiffusionBaseProblem, DomainSizeZ, 1.0);

SET_INT_PROP(DiffusionBaseProblem, CellsX, 250);
SET_INT_PROP(DiffusionBaseProblem, CellsY, 1);
SET_INT_PROP(DiffusionBaseProblem, CellsZ, 1);

// The default for the end time of the simulation
SET_SCALAR_PROP(DiffusionBaseProblem, EndTime, 1e6);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(DiffusionBaseProblem, InitialTimeStepSize, 1000);
}
} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 * \brief 1D problem which is driven by molecular diffusion.
 *
 * The domain is one meter long and completely filled with gas and
 * closed on all boundaries. Its left half exhibits a slightly higher
 * water concentration than the right one. After a while, the
 * concentration of water will be equilibrate due to molecular
 * diffusion.
 */
template <class TypeTag>
class DiffusionProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    DiffusionProblem(TimeManager &timeManager)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafGridView())
#else
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
#endif
    {
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
     * \copydoc VcfvProblem::name
     */
    static std::string name()
    { return std::string("diffusion_") + Model::name(); }

    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    { return K_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.35; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    { return materialParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     *
     * This problem sets no-flow boundaries everywhere.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    { values.setNoFlow(); }

    //! \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        if (onLeftSide_(pos))
            values.assignNaive(leftInitialFluidState_);
        else
            values.assignNaive(rightInitialFluidState_);
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
    bool onLeftSide_(const GlobalPosition &pos) const
    { return pos[0] < (this->boundingBoxMin()[0] + this->boundingBoxMax()[0]) / 2; }

    void setupInitialFluidStates_()
    {
        // create the initial fluid state for the left half of the domain
        leftInitialFluidState_.setTemperature(temperature_);

        Scalar Sl = 0.0;
        leftInitialFluidState_.setSaturation(lPhaseIdx, Sl);
        leftInitialFluidState_.setSaturation(gPhaseIdx, 1 - Sl);

        Scalar p = 1e5;
        leftInitialFluidState_.setPressure(lPhaseIdx, p);
        leftInitialFluidState_.setPressure(gPhaseIdx, p);

        Scalar xH2O = 0.01;
        leftInitialFluidState_.setMoleFraction(gPhaseIdx, H2OIdx, xH2O);
        leftInitialFluidState_.setMoleFraction(gPhaseIdx, N2Idx, 1 - xH2O);

        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        typename FluidSystem::ParameterCache paramCache;
        CFRP::solve(leftInitialFluidState_, paramCache, gPhaseIdx,
                    /*setViscosity=*/false, /*setEnthalpy=*/false);

        // create the initial fluid state for the right half of the domain
        rightInitialFluidState_.assign(leftInitialFluidState_);
        xH2O = 0.0;
        rightInitialFluidState_.setMoleFraction(gPhaseIdx, H2OIdx, xH2O);
        rightInitialFluidState_.setMoleFraction(gPhaseIdx, N2Idx, 1 - xH2O);
        CFRP::solve(rightInitialFluidState_, paramCache, gPhaseIdx,
                    /*setViscosity=*/false, /*setEnthalpy=*/false);
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> leftInitialFluidState_;
    Opm::CompositionalFluidState<Scalar, FluidSystem> rightInitialFluidState_;
    Scalar temperature_;
};

} // namespace Ewoms

#endif
