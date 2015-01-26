/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2010 by Melanie Darcis

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
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

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
                                 /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                 /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                 /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    typedef Opm::LinearMaterial<Traits> type;
};

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(ReservoirBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(ReservoirBaseProblem, EnableGravity, true);

// Reuse the last linearization if possible?
SET_BOOL_PROP(ReservoirBaseProblem, EnableLinearizationRecycling, true);

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
SET_STRING_PROP(ReservoirBaseProblem, GridFile, "data/reservoir.dgf");
}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup TestProblems
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
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) BlackOilFluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ReservoirProblem(Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);

        FluidSystem::initBegin();std::vector<std::pair<Scalar, Scalar> > Bo = {
            { 101353, 1.062 },
            { 1.82504e+06, 1.15 },
            { 3.54873e+06, 1.207 },
            { 6.99611e+06, 1.295 },
            { 1.38909e+07, 1.435 },
            { 1.73382e+07, 1.5 },
            { 2.07856e+07, 1.565 },
            { 2.76804e+07, 1.695 },
            { 3.45751e+07, 1.827 }
        };
        std::vector<std::pair<Scalar, Scalar> > muo = {
            { 101353, 0.00104 },
            { 1.82504e+06, 0.000975 },
            { 3.54873e+06, 0.00091 },
            { 6.99611e+06, 0.00083 },
            { 1.38909e+07, 0.000695 },
            { 1.73382e+07, 0.000641 },
            { 2.07856e+07, 0.000594 },
            { 2.76804e+07, 0.00051 },
            { 3.45751e+07, 0.000449 }
        };
        std::vector<std::pair<Scalar, Scalar> > Rs = {
            { 101353, 0.178108 },
            { 1.82504e+06, 16.1187 },
            { 3.54873e+06, 32.0594 },
            { 6.99611e+06, 66.0779 },
            { 1.38909e+07, 113.276 },
            { 1.73382e+07, 138.033 },
            { 2.07856e+07, 165.64 },
            { 2.76804e+07, 226.197 },
            { 3.45751e+07, 288.178 }
        };
        std::vector<std::pair<Scalar, Scalar> > Bg = {
            { 101353, 0.93576 },
            { 1.82504e+06, 0.0678972 },
            { 3.54873e+06, 0.0352259 },
            { 6.99611e+06, 0.0179498 },
            { 1.38909e+07, 0.00906194 },
            { 1.73382e+07, 0.00726527 },
            { 2.07856e+07, 0.00606375 },
            { 2.76804e+07, 0.00455343 },
            { 3.45751e+07, 0.00364386 },
            { 6.21542e+07, 0.00216723 }
        };
        std::vector<std::pair<Scalar, Scalar> > mug = {
            { 101353, 8e-06 },
            { 1.82504e+06, 9.6e-06 },
            { 3.54873e+06, 1.12e-05 },
            { 6.99611e+06, 1.4e-05 },
            { 1.38909e+07, 1.89e-05 },
            { 1.73382e+07, 2.08e-05 },
            { 2.07856e+07, 2.28e-05 },
            { 2.76804e+07, 2.68e-05 },
            { 3.45751e+07, 3.09e-05 },
            { 6.21542e+07, 4.7e-05 }
        };

        FluidSystem::setReferenceDensities(/*oil=*/786, /*water=*/1037, /*gas=*/0.97);
        FluidSystem::setGasFormationVolumeFactor(Bg);
        FluidSystem::setSaturatedOilGasDissolutionFactor(Rs);
        FluidSystem::setSaturatedOilFormationVolumeFactor(Bo);
        FluidSystem::setSaturatedOilViscosity(muo);
        FluidSystem::setGasViscosity(mug);
        FluidSystem::setWaterViscosity(9.6e-4);
        FluidSystem::setWaterCompressibility(1.450377e-10);
        FluidSystem::initEnd();

        pReservoir_ = 330e5;
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
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
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
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return EWOMS_GET_PARAM(TypeTag, std::string, SimulationName); }

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
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
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
     * \copydoc FvBaseMultiPhaseProblem::porosity
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
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
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
     * \copydoc FvBaseMultiPhaseProblem::temperature
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
     * \copydoc FvBaseProblem::boundary
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
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    { values.assignNaive(initialFluidState_); }

    /*!
     * \copydoc FvBaseProblem::constraints
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
        Scalar x = pos[0] - this->boundingBoxMin()[0];
        Scalar y = pos[dim - 1] - this->boundingBoxMin()[dim - 1];
        Scalar height = this->boundingBoxMax()[dim - 1] - this->boundingBoxMin()[dim - 1];
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        if ((onLeftBoundary_(pos) || onRightBoundary_(pos)) && y < height / 2) {
            // injectors
            auto fs = initialFluidState_;

            Scalar pInj = pReservoir_ * 1.5;
            fs.setPressure(waterPhaseIdx, pInj);
            fs.setPressure(oilPhaseIdx, pInj);
            fs.setPressure(gasPhaseIdx, pInj);
            fs.setSaturation(waterPhaseIdx, 1.0);
            fs.setSaturation(oilPhaseIdx, 0.0);
            fs.setSaturation(gasPhaseIdx, 0.0);

            // set the compositions to only water
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    fs.setMoleFraction(phaseIdx, compIdx, 0.0);

            // set the composition of the oil phase to the initial
            // composition
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(oilPhaseIdx, compIdx,
                                   initialFluidState_.moleFraction(oilPhaseIdx,
                                                                   compIdx));

            fs.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);

            constraints.setAllConstraint();
            constraints.assignNaive(fs);
        }
        else if (width / 2 - 1 < x && x < width / 2 + 1 && y > height / 2) {
            // producer
            auto fs = initialFluidState_;

            Scalar pProd = pReservoir_ / 1.5;
            fs.setPressure(waterPhaseIdx, pProd);
            fs.setPressure(oilPhaseIdx, pProd);
            fs.setPressure(gasPhaseIdx, pProd);
            fs.setSaturation(waterPhaseIdx, 0.0);
            fs.setSaturation(oilPhaseIdx, 1.0);
            fs.setSaturation(gasPhaseIdx, 0.0);

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
     * \copydoc FvBaseProblem::source
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
        fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar pw = pReservoir_;

        PhaseVector pC;
        const auto &matParams = fineMaterialParams_;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(oilPhaseIdx, pw + (pC[oilPhaseIdx] - pC[waterPhaseIdx]));
        fs.setPressure(waterPhaseIdx, pw + (pC[waterPhaseIdx] - pC[waterPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pw + (pC[gasPhaseIdx] - pC[waterPhaseIdx]));

        // reset all mole fractions to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(phaseIdx, compIdx, 0.0);

        //////
        // set composition of the gas and water phases
        //////
        fs.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
        fs.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);

        //////
        // set composition of the oil phase
        //////

        Scalar xoG = 0.95*FluidSystem::saturatedOilGasMoleFraction(fs.pressure(oilPhaseIdx));
        Scalar xoO = 1 - xoG;

        // finally set the oil-phase composition
        fs.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
        fs.setMoleFraction(oilPhaseIdx, oilCompIdx, xoO);
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

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
};
} // namespace Ewoms

#endif
