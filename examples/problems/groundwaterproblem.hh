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

#include <ewoms/models/immiscible/immiscibleproperties.hh>
#include <opm/simulators/linalg/parallelistlbackend.hh>

#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class GroundWaterProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(GroundWaterBaseProblem);

NEW_PROP_TAG(LensLowerLeftX);
NEW_PROP_TAG(LensLowerLeftY);
NEW_PROP_TAG(LensLowerLeftZ);
NEW_PROP_TAG(LensUpperRightX);
NEW_PROP_TAG(LensUpperRightY);
NEW_PROP_TAG(LensUpperRightZ);
NEW_PROP_TAG(Permeability);
NEW_PROP_TAG(PermeabilityLens);

SET_PROP(GroundWaterBaseProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(GroundWaterBaseProblem, Grid, Dune::YaspGrid<2>);
// SET_TYPE_PROP(GroundWaterBaseProblem, Grid, Dune::SGrid<2, 2>);

SET_TYPE_PROP(GroundWaterBaseProblem, Problem,
              Opm::GroundWaterProblem<TypeTag>);

SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftX, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftY, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftZ, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightX, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightY, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightZ, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, Permeability, 1e-10);
SET_SCALAR_PROP(GroundWaterBaseProblem, PermeabilityLens, 1e-12);

// Enable gravity
SET_BOOL_PROP(GroundWaterBaseProblem, EnableGravity, true);

// The default for the end time of the simulation
SET_SCALAR_PROP(GroundWaterBaseProblem, EndTime, 1);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(GroundWaterBaseProblem, InitialTimeStepSize, 1);

// The default DGF file to load
SET_STRING_PROP(GroundWaterBaseProblem, GridFile, "./data/groundwater_2d.dgf");

// Use the conjugated gradient linear solver with the default preconditioner (i.e.,
// ILU-0) from dune-istl
SET_TAG_PROP(GroundWaterBaseProblem, LinearSolverSplice, ParallelIstlLinearSolver);
SET_TYPE_PROP(GroundWaterBaseProblem, LinearSolverWrapper,
              Opm::Linear::SolverWrapperConjugatedGradients<TypeTag>);

END_PROPERTIES

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
class GroundWaterProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressure0Idx = Indices::pressure0Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

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

        lensLowerLeft_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftX);
        if (dim > 1)
            lensLowerLeft_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftY);
        if (dim > 2)
            lensLowerLeft_[2] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftY);

        lensUpperRight_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightX);
        if (dim > 1)
            lensUpperRight_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightY);
        if (dim > 2)
            lensUpperRight_[2] = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightY);

        intrinsicPerm_ = this->toDimMatrix_(EWOMS_GET_PARAM(TypeTag, Scalar, Permeability));
        intrinsicPermLens_ = this->toDimMatrix_(EWOMS_GET_PARAM(TypeTag, Scalar, PermeabilityLens));
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftX,
                             "The x-coordinate of the lens' lower-left corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightX,
                             "The x-coordinate of the lens' upper-right corner "
                             "[m].");

        if (dimWorld > 1) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftY,
                                 "The y-coordinate of the lens' lower-left "
                                 "corner [m].");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightY,
                                 "The y-coordinate of the lens' upper-right "
                                 "corner [m].");
        }

        if (dimWorld > 2) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftZ,
                                 "The z-coordinate of the lens' lower-left "
                                 "corner [m].");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightZ,
                                 "The z-coordinate of the lens' upper-right "
                                 "corner [m].");
        }

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Permeability,
                             "The intrinsic permeability [m^2] of the ambient "
                             "material.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, PermeabilityLens,
                             "The intrinsic permeability [m^2] of the lens.");
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
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return 273.15 + 10; } // 10C

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
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
                 const Context& context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        // const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);
        values[pressure0Idx] = 1.0e+5; // + 9.81*1.23*(20-globalPos[dim-1]);
    }

    /*!
     * \copydoc FvBaseProblem::source
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
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
