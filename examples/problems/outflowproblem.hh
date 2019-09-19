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
 * \copydoc Opm::OutflowProblem
 */
#ifndef EWOMS_OUTFLOW_PROBLEM_HH
#define EWOMS_OUTFLOW_PROBLEM_HH

#include <opm/models/pvs/pvsproperties.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/H2ON2LiquidPhaseFluidSystem.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
template <class TypeTag>
class OutflowProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(OutflowBaseProblem);

// Set the grid type
SET_TYPE_PROP(OutflowBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(OutflowBaseProblem, Problem, Opm::OutflowProblem<TypeTag>);

// Set fluid system
SET_PROP(OutflowBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    // Two-component single phase fluid system
    typedef Opm::H2ON2LiquidPhaseFluidSystem<Scalar> type;
};

// Disable gravity
SET_BOOL_PROP(OutflowBaseProblem, EnableGravity, false);

// Also write mass fractions to the output
SET_BOOL_PROP(OutflowBaseProblem, VtkWriteMassFractions, true);

// The default for the end time of the simulation
SET_SCALAR_PROP(OutflowBaseProblem, EndTime, 100);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(OutflowBaseProblem, InitialTimeStepSize, 1);

// The default DGF file to load
SET_STRING_PROP(OutflowBaseProblem, GridFile, "./data/outflow.dgf");

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where dissolved nitrogen is transported with the water
 *        phase from the left side to the right.
 *
 * The model domain is 1m times 1m and exhibits homogeneous soil
 * properties (\f$ \mathrm{K=10e-10, \Phi=0.4}\f$).  Initially the
 * domain is fully saturated by water without any nitrogen dissolved.
 *
 * At the left side, a free-flow condition defines a nitrogen mole
 * fraction of 0.02%.  The water phase flows from the left side to the
 * right due to the imposed pressure gradient of \f$1e5\,Pa/m\f$. The
 * nitrogen is transported with the water flow and leaves the domain
 * at the right boundary where an outflow boundary condition is
 * used.
 */
template <class TypeTag>
class OutflowProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // copy some indices for convenience
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = FluidSystem::numPhases,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    OutflowProblem(Simulator& simulator)
        : ParentType(simulator)
        , eps_(1e-6)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        temperature_ = 273.15 + 20;
        FluidSystem::init(/*minT=*/temperature_ - 1, /*maxT=*/temperature_ + 2,
                          /*numT=*/3,
                          /*minp=*/0.8e5, /*maxp=*/2.5e5, /*nump=*/500);

        // set parameters of porous medium
        perm_ = this->toDimMatrix_(1e-10);
        porosity_ = 0.4;
        tortuosity_ = 0.28;
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return "outflow"; }

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
     *
     * This problem assumes a temperature.
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; } // in [K]

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     *
     * This problem uses a constant intrinsic permeability.
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    { return perm_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     *
     * This problem uses a constant porosity.
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return porosity_; }

#if 0
    /*!
     * \brief Define the tortuosity \f$[?]\f$.
     *
     */
    template <class Context>
    Scalar tortuosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return tortuosity_; }

    /*!
     * \brief Define the dispersivity \f$[?]\f$.
     *
     */
    template <class Context>
    Scalar dispersivity(const Context& context,
                        unsigned spaceIdx, unsigned timeIdx) const
    { return 0; }
#endif

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

        if (onLeftBoundary_(globalPos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem,
                                         /*storeEnthalpy=*/false> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.setPressure(/*phaseIdx=*/0, fs.pressure(/*phaseIdx=*/0) + 1e5);

            Scalar xlN2 = 2e-4;
            fs.setMoleFraction(/*phaseIdx=*/0, N2Idx, xlN2);
            fs.setMoleFraction(/*phaseIdx=*/0, H2OIdx, 1 - xlN2);

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fs);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
                fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
            }

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onRightBoundary_(globalPos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem,
                                         /*storeEnthalpy=*/false> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);

            // impose an outflow boundary condition
            values.setOutFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            // no flow on top and bottom
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
        Opm::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        values.assignNaive(fs);
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

    template <class FluidState, class Context>
    void initialFluidState_(FluidState& fs, const Context& context,
                            unsigned spaceIdx, unsigned timeIdx) const
    {
        Scalar T = temperature(context, spaceIdx, timeIdx);
        // Scalar rho = FluidSystem::H2O::liquidDensity(T, /*pressure=*/1.5e5);
        // Scalar z = context.pos(spaceIdx, timeIdx)[dim - 1] -
        // this->boundingBoxMax()[dim - 1];
        // Scalar z = context.pos(spaceIdx, timeIdx)[dim - 1] -
        // this->boundingBoxMax()[dim - 1];

        fs.setSaturation(/*phaseIdx=*/0, 1.0);
        fs.setPressure(/*phaseIdx=*/0, 1e5 /* + rho*z */);
        fs.setMoleFraction(/*phaseIdx=*/0, H2OIdx, 1.0);
        fs.setMoleFraction(/*phaseIdx=*/0, N2Idx, 0);
        fs.setTemperature(T);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
            fs.setViscosity(phaseIdx, FluidSystem::viscosity(fs, paramCache, phaseIdx));
        }
    }

    const Scalar eps_;

    MaterialLawParams materialParams_;
    DimMatrix perm_;
    Scalar temperature_;
    Scalar porosity_;
    Scalar tortuosity_;
};
} // namespace Opm

#endif
