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
 * \copydoc Opm::ReservoirProblem
 */
#ifndef EWOMS_RESERVOIR_PROBLEM_HH
#define EWOMS_RESERVOIR_PROBLEM_HH

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <vector>
#include <string>

namespace Opm {
template <class TypeTag>
class ReservoirProblem;

} // namespace Opm

BEGIN_PROPERTIES


NEW_TYPE_TAG(ReservoirBaseProblem);

// Maximum depth of the reservoir
NEW_PROP_TAG(MaxDepth);
// The temperature inside the reservoir
NEW_PROP_TAG(Temperature);
// The width of producer/injector wells as a fraction of the width of the spatial domain
NEW_PROP_TAG(WellWidth);

// Set the grid type
SET_TYPE_PROP(ReservoirBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ReservoirBaseProblem, Problem, Opm::ReservoirProblem<TypeTag>);

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

// Enable constraint DOFs?
SET_BOOL_PROP(ReservoirBaseProblem, EnableConstraints, true);

// set the defaults for some problem specific properties
SET_SCALAR_PROP(ReservoirBaseProblem, MaxDepth, 2500);
SET_SCALAR_PROP(ReservoirBaseProblem, Temperature, 293.15);

//! The default for the end time of the simulation [s].
//!
//! By default this problem spans 1000 days (100 "settle down" days and 900 days of
//! production)
SET_SCALAR_PROP(ReservoirBaseProblem, EndTime, 1000.0*24*60*60);

// The default for the initial time step size of the simulation [s]
SET_SCALAR_PROP(ReservoirBaseProblem, InitialTimeStepSize, 100e3);

// The width of producer/injector wells as a fraction of the width of the spatial domain
SET_SCALAR_PROP(ReservoirBaseProblem, WellWidth, 0.01);

/*!
 * \brief Explicitly set the fluid system to the black-oil fluid system
 *
 * If the black oil model is used, this is superfluous because that model already sets
 * the FluidSystem property. Setting it explictly for the problem is a good idea anyway,
 * though because other models are more generic and thus do not assume a particular fluid
 * system.
 */
SET_PROP(ReservoirBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::BlackOilFluidSystem<Scalar> type;
};

// The default DGF file to load
SET_STRING_PROP(ReservoirBaseProblem, GridFile, "data/reservoir.dgf");

// increase the tolerance for this problem to get larger time steps
SET_SCALAR_PROP(ReservoirBaseProblem, NewtonTolerance, 1e-6);

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Some simple test problem for the black-oil VCVF discretization
 *        inspired by an oil reservoir.
 *
 * The domain is two-dimensional and exhibits a size of 6000m times 60m. Initially, the
 * reservoir is assumed by oil with a bubble point pressure of 20 MPa, which also the
 * initial pressure in the domain. No-flow boundaries are used for all boundaries. The
 * permeability of the lower 10 m is reduced compared to the upper 10 m of the domain
 * witch capillary pressure always being neglected. Three wells are approximated using
 * constraints: Two water-injector wells, one at the lower-left boundary one at the
 * lower-right boundary and one producer well in the upper part of the center of the
 * domain. The pressure for the producer is assumed to be 2/3 of the reservoir pressure,
 * the injector wells use a pressure which is 50% above the reservoir pressure.
 */
template <class TypeTag>
class ReservoirProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
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

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    typedef Opm::CompositionalFluidState<Scalar,
                                         FluidSystem,
                                         /*enableEnthalpy=*/true> InitialFluidState;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ReservoirProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        wellWidth_ = EWOMS_GET_PARAM(TypeTag, Scalar, WellWidth);

        std::vector<std::pair<Scalar, Scalar> > Bo = {
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

        Scalar rhoRefO = 786.0; // [kg]
        Scalar rhoRefG = 0.97; // [kg]
        Scalar rhoRefW = 1037.0; // [kg]
        FluidSystem::initBegin(/*numPvtRegions=*/1);
        FluidSystem::setEnableDissolvedGas(true);
        FluidSystem::setEnableVaporizedOil(false);
        FluidSystem::setReferenceDensities(rhoRefO, rhoRefW, rhoRefG, /*regionIdx=*/0);

        Opm::GasPvtMultiplexer<Scalar> *gasPvt = new Opm::GasPvtMultiplexer<Scalar>;
        gasPvt->setApproach(Opm::GasPvtMultiplexer<Scalar>::DryGasPvt);
        auto& dryGasPvt = gasPvt->template getRealPvt<Opm::GasPvtMultiplexer<Scalar>::DryGasPvt>();
        dryGasPvt.setNumRegions(/*numPvtRegion=*/1);
        dryGasPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
        dryGasPvt.setGasFormationVolumeFactor(/*regionIdx=*/0, Bg);
        dryGasPvt.setGasViscosity(/*regionIdx=*/0, mug);

        Opm::OilPvtMultiplexer<Scalar> *oilPvt = new Opm::OilPvtMultiplexer<Scalar>;
        oilPvt->setApproach(Opm::OilPvtMultiplexer<Scalar>::LiveOilPvt);
        auto& liveOilPvt = oilPvt->template getRealPvt<Opm::OilPvtMultiplexer<Scalar>::LiveOilPvt>();
        liveOilPvt.setNumRegions(/*numPvtRegion=*/1);
        liveOilPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
        liveOilPvt.setSaturatedOilGasDissolutionFactor(/*regionIdx=*/0, Rs);
        liveOilPvt.setSaturatedOilFormationVolumeFactor(/*regionIdx=*/0, Bo);
        liveOilPvt.setSaturatedOilViscosity(/*regionIdx=*/0, muo);

        Opm::WaterPvtMultiplexer<Scalar> *waterPvt = new Opm::WaterPvtMultiplexer<Scalar>;
        waterPvt->setApproach(Opm::WaterPvtMultiplexer<Scalar>::ConstantCompressibilityWaterPvt);
        auto& ccWaterPvt = waterPvt->template getRealPvt<Opm::WaterPvtMultiplexer<Scalar>::ConstantCompressibilityWaterPvt>();
        ccWaterPvt.setNumRegions(/*numPvtRegions=*/1);
        ccWaterPvt.setReferenceDensities(/*regionIdx=*/0, rhoRefO, rhoRefG, rhoRefW);
        ccWaterPvt.setViscosity(/*regionIdx=*/0, 9.6e-4);
        ccWaterPvt.setCompressibility(/*regionIdx=*/0, 1.450377e-10);

        gasPvt->initEnd();
        oilPvt->initEnd();
        waterPvt->initEnd();

        typedef std::shared_ptr<Opm::GasPvtMultiplexer<Scalar> > GasPvtSharedPtr;
        FluidSystem::setGasPvt(GasPvtSharedPtr(gasPvt));

        typedef std::shared_ptr<Opm::OilPvtMultiplexer<Scalar> > OilPvtSharedPtr;
        FluidSystem::setOilPvt(OilPvtSharedPtr(oilPvt));

        typedef std::shared_ptr<Opm::WaterPvtMultiplexer<Scalar> > WaterPvtSharedPtr;
        FluidSystem::setWaterPvt(WaterPvtSharedPtr(waterPvt));

        FluidSystem::initEnd();

        pReservoir_ = 330e5;
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-12);
        coarseK_ = this->toDimMatrix_(1e-11);

        // porosities
        finePorosity_ = 0.2;
        coarsePorosity_ = 0.3;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fineMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            fineMaterialParams_.setPcMaxSat(phaseIdx, 0.0);

            coarseMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            coarseMaterialParams_.setPcMaxSat(phaseIdx, 0.0);
        }

        // wrap up the initialization of the material law's parameters
        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        materialParams_.resize(this->model().numGridDof());
        ElementContext elemCtx(this->simulator());
        auto eIt = this->simulator().gridView().template begin<0>();
        const auto& eEndIt = this->simulator().gridView().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            elemCtx.updateStencil(*eIt);
            size_t nDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
            for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                const GlobalPosition& pos = elemCtx.pos(dofIdx, /*timeIdx=*/0);

                if (isFineMaterial_(pos))
                    materialParams_[globalDofIdx] = &fineMaterialParams_;
                else
                    materialParams_[globalDofIdx] = &coarseMaterialParams_;
            }
        }

        initFluidState_();

        // start the first ("settle down") episode for 100 days
        this->simulator().startNextEpisode(100.0*24*60*60);
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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, WellWidth,
                             "The width of producer/injector wells as a fraction of the width"
                             " of the spatial domain");
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return std::string("reservoir_") + Model::name() + "_" + Model::discretizationName(); }

    /*!
     * \copydoc FvBaseProblem::endEpisode
     */
    void endEpisode()
    {
        // in the second episode, the actual work is done (the first is "settle down"
        // episode). we need to use a pretty short initial time step here as the change
        // in conditions is quite abrupt.
        this->simulator().startNextEpisode(1e100);
        this->simulator().setTimeStepSize(5.0);
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
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     *
     * For this problem, a layer with high permability is located
     * above one with low permeability.
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return *materialParams_[globalIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned globalIdx) const
    { return *materialParams_[globalIdx]; }

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
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
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
    void boundary(BoundaryRateVector& values,
                  const Context& context OPM_UNUSED,
                  unsigned spaceIdx OPM_UNUSED,
                  unsigned timeIdx OPM_UNUSED) const
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
    void initial(PrimaryVariables& values,
                 const Context& context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        values.assignNaive(initialFluidState_);

#ifndef NDEBUG
        for (unsigned pvIdx = 0; pvIdx < values.size(); ++ pvIdx)
            assert(std::isfinite(values[pvIdx]));
#endif
    }

    /*!
     * \copydoc FvBaseProblem::constraints
     *
     * The reservoir problem places two water-injection wells on the lower-left and
     * lower-right of the domain and a production well in the middle. The injection wells
     * are fully water saturated with a higher pressure, the producer is fully oil
     * saturated with a lower pressure than the remaining reservoir.
     */
    template <class Context>
    void constraints(Constraints& constraintValues,
                     const Context& context,
                     unsigned spaceIdx,
                     unsigned timeIdx) const
    {
        if (this->simulator().episodeIndex() == 1)
            return; // no constraints during the "settle down" episode

        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (isInjector_(pos)) {
            constraintValues.setActive(true);
            constraintValues.assignNaive(injectorFluidState_);
        }
        else if (isProducer_(pos)) {
            constraintValues.setActive(true);
            constraintValues.assignNaive(producerFluidState_);
        }
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

    //! \}

private:
    void initFluidState_()
    {
        auto& fs = initialFluidState_;

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
        const auto& matParams = fineMaterialParams_;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(oilPhaseIdx, pw + (pC[oilPhaseIdx] - pC[waterPhaseIdx]));
        fs.setPressure(waterPhaseIdx, pw + (pC[waterPhaseIdx] - pC[waterPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pw + (pC[gasPhaseIdx] - pC[waterPhaseIdx]));

        // reset all mole fractions to 0
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(phaseIdx, compIdx, 0.0);

        //////
        // set composition of the gas and water phases
        //////
        fs.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
        fs.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);

        //////
        // set composition of the oil phase
        //////
        Scalar RsSat =
            FluidSystem::saturatedDissolutionFactor(fs, oilPhaseIdx, /*pvtRegionIdx=*/0);
        Scalar XoGSat = FluidSystem::convertRsToXoG(RsSat, /*pvtRegionIdx=*/0);
        Scalar xoGSat = FluidSystem::convertXoGToxoG(XoGSat, /*pvtRegionIdx=*/0);
        Scalar xoG = 0.95*xoGSat;
        Scalar xoO = 1.0 - xoG;

        // finally set the oil-phase composition
        fs.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
        fs.setMoleFraction(oilPhaseIdx, oilCompIdx, xoO);

        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        CFRP::solve(fs,
                    paramCache,
                    /*refPhaseIdx=*/oilPhaseIdx,
                    /*setViscosities=*/false,
                    /*setEnthalpies=*/false);

        // set up the fluid state used for the injectors
        auto& injFs = injectorFluidState_;
        injFs = initialFluidState_;

        Scalar pInj = pReservoir_ * 1.5;
        injFs.setPressure(waterPhaseIdx, pInj);
        injFs.setPressure(oilPhaseIdx, pInj);
        injFs.setPressure(gasPhaseIdx, pInj);
        injFs.setSaturation(waterPhaseIdx, 1.0);
        injFs.setSaturation(oilPhaseIdx, 0.0);
        injFs.setSaturation(gasPhaseIdx, 0.0);

        // set the composition of the phases to immiscible
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                injFs.setMoleFraction(phaseIdx, compIdx, 0.0);

        injFs.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        injFs.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);
        injFs.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);

        CFRP::solve(injFs,
                    paramCache,
                    /*refPhaseIdx=*/waterPhaseIdx,
                    /*setViscosities=*/true,
                    /*setEnthalpies=*/false);

        // set up the fluid state used for the producer
        auto& prodFs = producerFluidState_;
        prodFs = initialFluidState_;

        Scalar pProd = pReservoir_ / 1.5;
        prodFs.setPressure(waterPhaseIdx, pProd);
        prodFs.setPressure(oilPhaseIdx, pProd);
        prodFs.setPressure(gasPhaseIdx, pProd);
        prodFs.setSaturation(waterPhaseIdx, 0.0);
        prodFs.setSaturation(oilPhaseIdx, 1.0);
        prodFs.setSaturation(gasPhaseIdx, 0.0);

        CFRP::solve(prodFs,
                    paramCache,
                    /*refPhaseIdx=*/oilPhaseIdx,
                    /*setViscosities=*/true,
                    /*setEnthalpies=*/false);
    }

    bool isProducer_(const GlobalPosition& pos) const
    {
        Scalar x = pos[0] - this->boundingBoxMin()[0];
        Scalar y = pos[dim - 1] - this->boundingBoxMin()[dim - 1];
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar height = this->boundingBoxMax()[dim - 1] - this->boundingBoxMin()[dim - 1];

        // only the upper half of the center section of the spatial domain is assumed to
        // be the producer
        if (y <= height/2.0)
            return false;

        return width/2.0 - width*1e-5 < x && x < width/2.0 + width*(wellWidth_ + 1e-5);
    }

    bool isInjector_(const GlobalPosition& pos) const
    {
        Scalar x = pos[0] - this->boundingBoxMin()[0];
        Scalar y = pos[dim - 1] - this->boundingBoxMin()[dim - 1];
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar height = this->boundingBoxMax()[dim - 1] - this->boundingBoxMin()[dim - 1];

        // only the lower half of the leftmost and rightmost part of the spatial domain
        // are assumed to be the water injectors
        if (y > height/2.0)
            return false;

        return x < width*wellWidth_ - width*1e-5 || x > width*(1.0 - wellWidth_) + width*1e-5;
    }

    bool isFineMaterial_(const GlobalPosition& pos) const
    { return pos[dim - 1] > layerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar layerBottom_;
    Scalar pReservoir_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;
    std::vector<const MaterialLawParams*> materialParams_;

    InitialFluidState initialFluidState_;
    InitialFluidState injectorFluidState_;
    InitialFluidState producerFluidState_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar wellWidth_;
};
} // namespace Opm

#endif
