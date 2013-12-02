// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \copydoc Ewoms::LensProblem
 */
#ifndef EWOMS_LENS_PROBLEM_HH
#define EWOMS_LENS_PROBLEM_HH

#include "lensgridcreator.hh"

#include <ewoms/models/immiscible/immiscibleproperties.hh>
#include <ewoms/linear/parallelamgbackend.hh>

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Dnapl.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <iostream>

namespace Ewoms {

template <class TypeTag>
class LensProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(LensBaseProblem);

// declare the properties specific for the lens problem
NEW_PROP_TAG(LensLowerLeftX);
NEW_PROP_TAG(LensLowerLeftY);
NEW_PROP_TAG(LensLowerLeftZ);
NEW_PROP_TAG(LensUpperRightX);
NEW_PROP_TAG(LensUpperRightY);
NEW_PROP_TAG(LensUpperRightZ);

// set the GridCreator property
SET_TYPE_PROP(LensBaseProblem, GridCreator, Ewoms::LensGridCreator<TypeTag>);

// Retrieve the grid type from the grid creator
SET_TYPE_PROP(LensBaseProblem, Grid,
              typename GET_PROP_TYPE(TypeTag, GridCreator)::Grid);

// Set the problem property
SET_TYPE_PROP(LensBaseProblem, Problem, Ewoms::LensProblem<TypeTag>);

// Set the wetting phase
SET_PROP(LensBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(LensBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::DNAPL<Scalar> > type;
};

// Set the material Law
SET_PROP(LensBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nPhaseIdx>
    Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedVanGenuchten<Traits> EffectiveLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffectiveLaw> type;
};

// Use the algebraic multi-grid linear solver for this problem
SET_TAG_PROP(LensBaseProblem, LinearSolver, ParallelAmgBackend);

// Enable partial reassembly of the jacobian matrix?
// SET_BOOL_PROP(LensBaseProblem, EnablePartialReassemble, true);

// Enable reuse of jacobian matrices?
// SET_BOOL_PROP(LensBaseProblem, EnableJacobianRecycling, true);

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(LensBaseProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(LensBaseProblem, NumericDifferenceMethod, +1);

// Enable gravity
SET_BOOL_PROP(LensBaseProblem, EnableGravity, true);

// define the properties specific for the lens problem
SET_SCALAR_PROP(LensBaseProblem, LensLowerLeftX, 1.0);
SET_SCALAR_PROP(LensBaseProblem, LensLowerLeftY, 2.0);
SET_SCALAR_PROP(LensBaseProblem, LensLowerLeftZ, 0.0);
SET_SCALAR_PROP(LensBaseProblem, LensUpperRightX, 4.0);
SET_SCALAR_PROP(LensBaseProblem, LensUpperRightY, 3.0);
SET_SCALAR_PROP(LensBaseProblem, LensUpperRightZ, 1.0);

SET_SCALAR_PROP(LensBaseProblem, DomainSizeX, 6.0);
SET_SCALAR_PROP(LensBaseProblem, DomainSizeY, 4.0);
SET_SCALAR_PROP(LensBaseProblem, DomainSizeZ, 1.0);

SET_INT_PROP(LensBaseProblem, CellsX, 48);
SET_INT_PROP(LensBaseProblem, CellsY, 32);
SET_INT_PROP(LensBaseProblem, CellsZ, 16);

// The default for the end time of the simulation
SET_SCALAR_PROP(LensBaseProblem, EndTime, 30e3);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(LensBaseProblem, InitialTimeStepSize, 250);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the model, the depth of the domain is implicitly
 * assumed to be 1 m everywhere.
 *
 * On the top and the bottom of the domain no-flow boundary conditions
 * are used, while free-flow conditions apply on the left and right
 * boundaries; DNAPL is injected at the top boundary from 3m to 4m at
 * a rate of 0.04 kg/(s m^2).
 *
 * At the boundary on the left, a free-flow condition using the
 * hydrostatic pressure scaled by a factor of 1.125 is imposed, while
 * on the right, it is just the hydrostatic pressure. The DNAPL
 * saturation on both sides is zero.
 */
template <class TypeTag>
class LensProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag,
                                   BoundaryRateVector) BoundaryRateVector;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    LensProblem(TimeManager &timeManager)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafGridView())
#else
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
#endif
    {
        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C
        lensLowerLeft_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftX);
        lensLowerLeft_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftY);
        lensUpperRight_[0] = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightX);
        lensUpperRight_[1] = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightY);

        if (dimWorld == 3) {
            lensLowerLeft_[2] = EWOMS_GET_PARAM(TypeTag, Scalar, LensLowerLeftZ);
            lensUpperRight_[2]
                = EWOMS_GET_PARAM(TypeTag, Scalar, LensUpperRightZ);
        }

        // parameters for the Van Genuchten law
        // alpha and n
        lensMaterialParams_.setVgAlpha(0.00045);
        lensMaterialParams_.setVgN(7.3);
        outerMaterialParams_.setVgAlpha(0.0037);
        outerMaterialParams_.setVgN(4.7);

        lensMaterialParams_.finalize();
        outerMaterialParams_.finalize();

        lensK_ = this->toDimMatrix_(9.05e-12);
        outerK_ = this->toDimMatrix_(4.6e-10);

        if (dimWorld == 3) {
            this->gravity_ = 0;
            this->gravity_[1] = -9.81;
        }
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftX,
                             "The x-coordinate of the lens' lower-left corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftY,
                             "The y-coordinate of the lens' lower-left corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightX,
                             "The x-coordinate of the lens' upper-right corner "
                             "[m].");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightY,
                             "The y-coordinate of the lens' upper-right corner "
                             "[m].");

        if (dimWorld == 3) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensLowerLeftZ,
                                 "The z-coordinate of the lens' lower-left "
                                 "corner [m].");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, LensUpperRightZ,
                                 "The z-coordinate of the lens' upper-right "
                                 "corner [m].");
        }
    };

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.4; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    //! \}

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "lens_" << this->model().name();
        return oss.str();
    }

    /*!
     * \copydoc VcfvProblem::postTimeStep
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: " << storage << std::endl;
        }
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            // free flow boundary
            Scalar densityW
                = WettingPhase::density(temperature_, /*pressure=*/1e5);

            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;

            // set wetting phase pressure and saturation
            if (onLeftBoundary_(pos)) {
                Scalar height = this->bboxMax()[1] - this->bboxMin()[1];
                Scalar depth = this->bboxMax()[1] - pos[1];
                Scalar alpha = (1 + 1.5 / height);

                // hydrostatic pressure scaled by alpha
                pw = 1e5 - alpha * densityW * this->gravity()[1] * depth;
                Sw = 1.0;
            }
            else {
                Scalar depth = this->bboxMax()[1] - pos[1];

                // hydrostatic pressure
                pw = 1e5 - densityW * this->gravity()[1] * depth;
                Sw = 1.0;
            }

            // specify a full fluid state using pw and Sw
            const MaterialLawParams &matParams
                = this->materialLawParams(context, spaceIdx, timeIdx);

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                                      /*storeEnthalpy=*/false> fs;
            fs.setSaturation(wPhaseIdx, Sw);
            fs.setSaturation(nPhaseIdx, 1 - Sw);
            fs.setTemperature(T);

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wPhaseIdx, pw);
            fs.setPressure(nPhaseIdx, pw + pC[nPhaseIdx] - pC[wPhaseIdx]);

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiNEqIdx] = -0.04; // kg / (m^2 * s)

            // impose a forced flow boundary
            values.setMassRate(massRate);
        }
        else {
            // no flow boundary
            values.setNoFlow();
        }
    }

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
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        Scalar depth = this->bboxMax()[1] - pos[1];

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = 1.0;
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fs, wPhaseIdx);
        Scalar densityW = FluidSystem::density(fs, paramCache, wPhaseIdx);

        // hydrostatic pressure (assuming incompressibility)
        Scalar pw = 1e5 - densityW * this->gravity()[1] * depth;

        // calculate the capillary pressure
        const MaterialLawParams &matParams
            = this->materialLawParams(context, spaceIdx, timeIdx);
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // make a full fluid state
        fs.setPressure(wPhaseIdx, pw);
        fs.setPressure(nPhaseIdx, pw + (pC[wPhaseIdx] - pC[nPhaseIdx]));

        // assign the primary variables
        values.assignNaive(fs);
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
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] - eps_ || pos[i] > lensUpperRight_[i]
                                                              + eps_)
                return false;
        }
        return true;
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - pos[0]) / width;
        return onUpperBoundary_(pos) && 0.5 < lambda && lambda < 2.0 / 3.0;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    DimMatrix lensK_;
    DimMatrix outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;

    Scalar temperature_;
    Scalar eps_;
};

} // namespace Ewoms

#endif
