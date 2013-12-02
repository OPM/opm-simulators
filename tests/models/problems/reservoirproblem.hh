// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \copydoc Ewoms::ReservoirProblem
 */
#ifndef EWOMS_RESERVOIR_PROBLEM_HH
#define EWOMS_RESERVOIR_PROBLEM_HH

#include <ewoms/models/blackoil/blackoilproperties.hh>

#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <vector>
#include <string>

namespace Ewoms {
template <class TypeTag>
class ReservoirProblem;
}

namespace Opm {
namespace Properties {

NEW_TYPE_TAG(ReservoirBaseProblem);

// Problem specific properties:

// Maximum depth of the reservoir
NEW_PROP_TAG(MaxDepth);
// The temperature inside the reservoir
NEW_PROP_TAG(Temperature);
// The name of the simulation (used for writing VTK files)
NEW_PROP_TAG(SimulationName);

// Set the grid type
SET_TYPE_PROP(ReservoirBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ReservoirBaseProblem, Problem, Ewoms::ReservoirProblem<TypeTag>);

// Set the material Law
SET_PROP(ReservoirBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Opm::
        ThreePhaseMaterialTraits<Scalar,
                                 /*wettingPhaseIdx=*/FluidSystem::wPhaseIdx,
                                 /*nonWettingPhaseIdx=*/FluidSystem::oPhaseIdx,
                                 /*gasPhaseIdx=*/FluidSystem::gPhaseIdx> Traits;

public:
    typedef Opm::LinearMaterial<Traits> type;
};

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(ReservoirBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(ReservoirBaseProblem, EnableGravity, true);

// Reuse Jacobian matrices if possible?
SET_BOOL_PROP(ReservoirBaseProblem, EnableJacobianRecycling, true);

// Smoothen the upwinding method?
SET_BOOL_PROP(ReservoirBaseProblem, EnableSmoothUpwinding, false);

// Enable constraint DOFs?
SET_BOOL_PROP(ReservoirBaseProblem, EnableConstraints, true);

// set the defaults for some problem specific properties
SET_SCALAR_PROP(ReservoirBaseProblem, MaxDepth, 2500);
SET_SCALAR_PROP(ReservoirBaseProblem, Temperature, 293.15);
SET_STRING_PROP(ReservoirBaseProblem, SimulationName, "reservoir");

// The default for the end time of the simulation [s]
SET_SCALAR_PROP(ReservoirBaseProblem, EndTime, 100);

// The default for the initial time step size of the simulation [s]
SET_SCALAR_PROP(ReservoirBaseProblem, InitialTimeStepSize, 10);

// The default DGF file to load
SET_STRING_PROP(ReservoirBaseProblem, GridFile, "grids/reservoir.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 *
 * \brief Some simple test problem for the black-oil VCVF discretization
 *        inspired by an oil reservoir.
 *
 * The domain is two-dimensional and exhibits a size of 6000m times
 * 60m. Initially, the reservoir is assumed by oil with a bubble point
 * pressure of 20 MPa, which also the initial pressure in the
 * domain. No-flow boundaries are used for all boundaries. The
 * permeability of the lower 10 m is reduced compared to the upper 10
 * m of the domain witch capillary pressure always being
 * neglected. Three wells are approximated using constraints: Two
 * water-injector wells, one at the lower-left boundary one at the
 * lower-right boundary and one producer well in the upper part of the
 * center of the domain. The pressure for the producer is assumed to
 * be 2/3 of the reservoir pressure, the injector wells use a pressure
 * which is 50% above the reservoir pressure.
 */
template <class TypeTag>
class ReservoirProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };
    enum { oPhaseIdx = FluidSystem::oPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { gCompIdx = FluidSystem::gCompIdx };
    enum { oCompIdx = FluidSystem::oCompIdx };
    enum { wCompIdx = FluidSystem::wCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag,
                                   BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag,
                                   BlackOilFluidState) BlackOilFluidState;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ReservoirProblem(TimeManager &timeManager)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafGridView())
#else
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
#endif
    {
        eps_ = 1e-6;

        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        name_ = EWOMS_GET_PARAM(TypeTag, std::string, SimulationName);

        FluidSystem::initBegin();
        std::vector<std::pair<Scalar, Scalar> > Bg
            = { { 1.013529e+05, 9.998450e-01 },
                { 2.757903e+06, 3.075500e-02 },
                { 5.515806e+06, 1.537947e-02 },
                { 8.273709e+06, 1.021742e-02 },
                { 1.103161e+07, 7.662783e-03 },
                { 1.378951e+07, 6.151899e-03 },
                { 1.654742e+07, 5.108709e-03 },
                { 1.930532e+07, 4.378814e-03 },
                { 2.206322e+07, 3.857780e-03 },
                { 2.482113e+07, 3.388401e-03 },
                { 2.757903e+07, 3.049842e-03 } };
        std::vector<std::pair<Scalar, Scalar> > Bo
            = { { 1.013529e+05, 1.000000e+00 },
                { 2.757903e+06, 1.012000e+00 },
                { 5.515806e+06, 1.025500e+00 },
                { 8.273709e+06, 1.038000e+00 },
                { 1.103161e+07, 1.051000e+00 },
                { 1.378951e+07, 1.063000e+00 },
                { 1.654742e+07, 1.075000e+00 },
                { 1.930532e+07, 1.087000e+00 },
                { 2.206322e+07, 1.098500e+00 },
                { 2.482113e+07, 1.110000e+00 },
                { 2.757903e+07, 1.120000e+00 } };
        std::vector<std::pair<Scalar, Scalar> > Rs
            = { { 1.013529e+05, 0.000000e+00 },
                { 2.757903e+06, 2.938776e+01 },
                { 5.515806e+06, 5.966605e+01 },
                { 8.273709e+06, 8.905380e+01 },
                { 1.103161e+07, 1.184416e+02 },
                { 1.378951e+07, 1.474731e+02 },
                { 1.654742e+07, 1.754360e+02 },
                { 1.930532e+07, 2.012616e+02 },
                { 2.206322e+07, 2.261967e+02 },
                { 2.482113e+07, 2.475696e+02 },
                { 2.757903e+07, 2.671614e+02 } };
        std::vector<std::pair<Scalar, Scalar> > muo
            = { { 1.013529e+05, 1.200000e-03 },
                { 2.757903e+06, 1.170000e-03 },
                { 5.515806e+06, 1.140000e-03 },
                { 8.273709e+06, 1.110000e-03 },
                { 1.103161e+07, 1.080000e-03 },
                { 1.378951e+07, 1.060000e-03 },
                { 1.654742e+07, 1.030000e-03 },
                { 1.930532e+07, 1.000000e-03 },
                { 2.206322e+07, 9.800000e-04 },
                { 2.482113e+07, 9.500000e-04 },
                { 2.757903e+07, 9.400000e-04 } };
        std::vector<std::pair<Scalar, Scalar> > mug
            = { { 1.013529e+05, 1.250000e-05 },
                { 2.757903e+06, 1.300000e-05 },
                { 5.515806e+06, 1.350000e-05 },
                { 8.273709e+06, 1.400000e-05 },
                { 1.103161e+07, 1.450000e-05 },
                { 1.378951e+07, 1.500000e-05 },
                { 1.654742e+07, 1.550000e-05 },
                { 1.930532e+07, 1.600000e-05 },
                { 2.206322e+07, 1.650000e-05 },
                { 2.482113e+07, 1.700000e-05 },
                { 2.757903e+07, 1.750000e-05 }, };
        FluidSystem::setGasFormationVolumeFactor(Bg);
        FluidSystem::setOilFormationVolumeFactor(Bo);
        FluidSystem::setGasDissolutionFactor(Rs);
        FluidSystem::setOilViscosity(muo);
        FluidSystem::setGasViscosity(mug);
        FluidSystem::setWaterViscosity(9.6e-4);
        FluidSystem::setWaterCompressibility(1.450377e-10);
        FluidSystem::setSurfaceDensities(/*oil=*/720.51,
                                         /*water=*/1009.32,
                                         /*gas=*/1.1245);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            FluidSystem::setReferenceVolumeFactor(phaseIdx, 1.0);
        FluidSystem::initEnd();

        pReservoir_ = 20e6;
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-12);
        coarseK_ = this->toDimMatrix_(1e-11);

        // porosities
        finePorosity_ = 0.2;
        coarsePorosity_ = 0.3;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fineMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            fineMaterialParams_.setPcMaxSat(phaseIdx, 0.0);

            coarseMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            coarseMaterialParams_.setPcMaxSat(phaseIdx, 0.0);
        }

        // wrap up the initialization of the material law's parameters
        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        initFluidState_();
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     *
     * For this problem, a layer with high permability is located
     * above one with low permeability.
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
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \copydoc VcfvMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    const std::string name() const
    { return name_; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     *
     * The black-oil model assumes constant temperature to define its
     * parameters. Although temperature is thus not really used by the
     * model, it gets written to the VTK output. Who nows, maybe we
     * will need it one day?
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     *
     * The reservoir problem uses constraints to approximate
     * extraction and production wells, so all boundaries are no-flow.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        // no flow on top and bottom
        values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volume terms
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        //////
        // set the primary variables
        //////
        values.assignNaive(initialFluidState_);
    }

    /*!
     * \copydoc VcfvProblem::constraints
     *
     * The reservoir problem places two water-injection wells on the
     * lower parts of the left and right edges of the domains and on
     * production well in the middle. The injection wells are fully
     * water saturated with a higher pressure, the producer is fully
     * oil saturated with a lower pressure than the remaining
     * reservoir.
     */
    template <class Context>
    void constraints(Constraints &constraints, const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        Scalar x = pos[0] - this->bboxMin()[0];
        Scalar y = pos[dim - 1] - this->bboxMin()[dim - 1];
        Scalar height = this->bboxMax()[dim - 1] - this->bboxMin()[dim - 1];
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        if ((onLeftBoundary_(pos) || onRightBoundary_(pos)) && y < height / 2) {
            // injectors
            auto fs = initialFluidState_;

            Scalar pInj = pReservoir_ * 1.5;
            fs.setPressure(wPhaseIdx, pInj);
            fs.setPressure(oPhaseIdx, pInj);
            fs.setPressure(gPhaseIdx, pInj);
            fs.setSaturation(wPhaseIdx, 1.0);
            fs.setSaturation(oPhaseIdx, 0.0);
            fs.setSaturation(gPhaseIdx, 0.0);

            // set the compositions to only water
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    fs.setMoleFraction(phaseIdx, compIdx, 0.0);

            // set the composition of the oil phase to the initial
            // composition
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(oPhaseIdx, compIdx,
                                   initialFluidState_.moleFraction(oPhaseIdx,
                                                                   compIdx));

            fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);

            constraints.setAllConstraint();
            constraints.assignNaive(fs);
        }
        else if (width / 2 - 1 < x && x < width / 2 + 1 && y > height / 2) {
            // producer
            auto fs = initialFluidState_;

            Scalar pProd = pReservoir_ / 1.5;
            fs.setPressure(wPhaseIdx, pProd);
            fs.setPressure(oPhaseIdx, pProd);
            fs.setPressure(gPhaseIdx, pProd);
            fs.setSaturation(wPhaseIdx, 0.0);
            fs.setSaturation(oPhaseIdx, 1.0);
            fs.setSaturation(gPhaseIdx, 0.0);

            // set the compositions to the initial composition
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    fs.setMoleFraction(phaseIdx, compIdx,
                                       initialFluidState_.moleFraction(phaseIdx,
                                                                       compIdx));

            constraints.setAllConstraint();
            constraints.assignNaive(fs);
        }
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
    void initFluidState_()
    {
        auto &fs = initialFluidState_;

        //////
        // set temperatures
        //////
        fs.setTemperature(temperature_);

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::oPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::wPhaseIdx, 0.0);
        fs.setSaturation(FluidSystem::gPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar pw = pReservoir_;

        PhaseVector pC;
        const auto &matParams = fineMaterialParams_;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(oPhaseIdx, pw + (pC[oPhaseIdx] - pC[wPhaseIdx]));
        fs.setPressure(wPhaseIdx, pw + (pC[wPhaseIdx] - pC[wPhaseIdx]));
        fs.setPressure(gPhaseIdx, pw + (pC[gPhaseIdx] - pC[wPhaseIdx]));

        // reset all mole fractions to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(phaseIdx, compIdx, 0.0);

        //////
        // set composition of the gas and water phases
        //////
        fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        fs.setMoleFraction(gPhaseIdx, gCompIdx, 1.0);

        //////
        // set composition of the oil phase
        //////

        // retrieve the relevant black-oil parameters from the fluid
        // system.
        Scalar pSat = pReservoir_; // the saturation pressure of the oil
        Scalar Bo = FluidSystem::oilFormationVolumeFactor(pSat);
        Scalar Rs = FluidSystem::gasDissolutionFactor(pSat);
        Scalar rhoo = FluidSystem::surfaceDensity(oPhaseIdx) / Bo;
        Scalar rhogref = FluidSystem::surfaceDensity(gPhaseIdx);

        // calculate composition of oil phase in terms of mass
        // fractions.
        Scalar XoG = Rs * rhogref / rhoo;

        // convert mass to mole fractions
        Scalar MG = FluidSystem::molarMass(gCompIdx);
        Scalar MO = FluidSystem::molarMass(oCompIdx);

        Scalar xoG = XoG * MO / ((MO - MG) * XoG + MG);
        Scalar xoO = 1 - xoG;

        // finally set the oil-phase composition
        fs.setMoleFraction(oPhaseIdx, gCompIdx, xoG);
        fs.setMoleFraction(oPhaseIdx, oCompIdx, xoO);
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onRightBoundary_(pos) && (5 < pos[1]) && (pos[1] < 15); }

    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim - 1] > layerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar layerBottom_;
    Scalar pReservoir_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    BlackOilFluidState initialFluidState_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    std::string name_;
};
} // namespace Ewoms

#endif
