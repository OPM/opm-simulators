// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Ewoms::InfiltrationProblem
 */
#ifndef EWOMS_INFILTRATION_PROBLEM_HH
#define EWOMS_INFILTRATION_PROBLEM_HH

#include <ewoms/models/pvs/pvsproperties.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/3p/3pParkerVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/3pAdapter.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/heatconduction/Somerton.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class InfiltrationProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(InfiltrationBaseProblem);

// Set the grid type
SET_TYPE_PROP(InfiltrationBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationBaseProblem, Problem,
              Ewoms::InfiltrationProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(
    InfiltrationBaseProblem, FluidSystem,
    Opm::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity?
SET_BOOL_PROP(InfiltrationBaseProblem, EnableGravity, true);

// Write newton convergence?
SET_BOOL_PROP(InfiltrationBaseProblem, NewtonWriteConvergence, false);

// Maximum tolerated relative error in the Newton method
SET_SCALAR_PROP(InfiltrationBaseProblem, NewtonRelativeTolerance, 1e-8);

// -1 backward differences, 0: central differences, +1: forward differences
SET_INT_PROP(InfiltrationBaseProblem, NumericDifferenceMethod, 1);

// Set the material Law
SET_PROP(InfiltrationBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    // define the three-phase material law
    typedef Opm::ThreePParkerVanGenuchten<Scalar> ThreePLaw;

public:
    // wrap the three-phase law in an adaptor to make use the generic
    // material law API
    typedef Opm::ThreePAdapter<wPhaseIdx, nPhaseIdx, gPhaseIdx, ThreePLaw> type;
};

// Set the heat conduction law
SET_PROP(InfiltrationBaseProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::Somerton<FluidSystem, Scalar> type;
};

// The default for the end time of the simulation
SET_SCALAR_PROP(InfiltrationBaseProblem, EndTime, 6e3);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(InfiltrationBaseProblem, InitialTimeStepSize, 60);

// The default DGF file to load
SET_STRING_PROP(InfiltrationBaseProblem, GridFile,
                "./grids/infiltration_50x3.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
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
class InfiltrationProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag,
                                   BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
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
        wPhaseIdx = FluidSystem::wPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    InfiltrationProblem(TimeManager &timeManager)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafGridView()),
#else
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView()),
#endif
          eps_(1e-6)
    {
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

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgN(4.);
        materialParams_.setkrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);

        materialParams_.checkDefined();
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::shouldWriteRestartFile
     *
     * This problem writes a restart file after every time step.
     */
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \copydoc VcfvProblem::name
     */
    const std::string name() const
    {
        std::ostringstream oss;
        oss << "infiltration_" << this->model().name();
        return oss.str();
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        // const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        // if (isFineMaterial_(pos))
        //     return finePorosity_;
        // else
        //     return coarsePorosity_;
        return porosity_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    { return materialParams_; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::heatCapacitySolid
     *
     * In this case, we assume the rock-matrix to be quartz.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx,
                             int timeIdx) const
    {
        return 850.     // specific heat capacity [J / (kg K)]
               * 2650.; // density of sand [kg/m^3]
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
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector molarRate(0.0);
            molarRate[conti0EqIdx + NAPLIdx] = -0.001;

            values.setMolarRate(molarRate);
            Valgrind::CheckDefined(values);
        }
        else
            values.setNoFlow();
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
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        Valgrind::CheckDefined(values);
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
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onUpperBoundary_(pos) && 50 < pos[0] && pos[0] < 75; }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState &fs, const Context &context,
                            int spaceIdx, int timeIdx) const
    {
        const GlobalPosition pos = context.pos(spaceIdx, timeIdx);
        Scalar y = pos[1];
        Scalar x = pos[0];

        Scalar densityW = 1000.0;
        Scalar pc = 9.81 * densityW * (y - (5 - 5e-4 * x));
        if (pc < 0.0)
            pc = 0.0;

        // set pressures
        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        Scalar Sw = invertPCGW_(pc, matParams);
        Scalar Swr = matParams.satResidual(wPhaseIdx);
        Scalar Sgr = matParams.satResidual(gPhaseIdx);
        if (Sw < Swr)
            Sw = Swr;
        if (Sw > 1 - Sgr)
            Sw = 1 - Sgr;
        Scalar Sg = 1 - Sw;

        Valgrind::CheckDefined(Sw);
        Valgrind::CheckDefined(Sg);

        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(gPhaseIdx, Sg);
        fs.setSaturation(nPhaseIdx, 0);

        // set temperature of all phases
        fs.setTemperature(temperature_);

        // compute pressures
        Scalar pcAll[numPhases];
        Scalar pg = 1e5;
        if (onLeftBoundary_(pos))
            pg += 10e3;
        MaterialLaw::capillaryPressures(pcAll, matParams, fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fs.setPressure(phaseIdx, pg + (pcAll[phaseIdx] - pcAll[gPhaseIdx]));

        // set composition of gas phase
        fs.setMoleFraction(gPhaseIdx, H2OIdx, 1e-6);
        fs.setMoleFraction(gPhaseIdx, airIdx,
                           1 - fs.moleFraction(gPhaseIdx, H2OIdx));
        fs.setMoleFraction(gPhaseIdx, NAPLIdx, 0);

        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        typename FluidSystem::ParameterCache paramCache;
        CFRP::solve(fs, paramCache, gPhaseIdx,
                    /*setViscosity=*/false,
                    /*setEnthalpy=*/false);

        fs.setMoleFraction(wPhaseIdx, H2OIdx,
                           1 - fs.moleFraction(wPhaseIdx, H2OIdx));
    }

    static Scalar invertPCGW_(Scalar pcIn, const MaterialLawParams &pcParams)
    {
        Scalar lower, upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar Sw, pcGW;
        lower = 0.0;
        upper = 1.0;
        for (k = 1; k <= 25; k++) {
            Sw = 0.5 * (upper + lower);
            pcGW = MaterialLaw::pCGW(pcParams, Sw);
            Scalar delta = pcGW - pcIn;
            if (delta < 0.)
                delta *= -1.;
            if (delta < bisLimit) {
                return (Sw);
            }
            if (k == maxIt) {
                return (Sw);
            }
            if (pcGW > pcIn)
                lower = Sw;
            else
                upper = Sw;
        }
        return 0;
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        return 70. <= pos[0] && pos[0] <= 85. && 7.0 <= pos[1] && pos[1] <= 7.50;
    }

    DimMatrix fineK_;
    DimMatrix coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Ewoms

#endif
