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
 * \copydoc Opm::FractureProblem
 */
#ifndef EWOMS_FRACTURE_PROBLEM_HH
#define EWOMS_FRACTURE_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
// avoid reordering of macro elements, otherwise this problem won't work
#define DISABLE_ALUGRID_SFC_ORDERING 1
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#else
#error "dune-alugrid not found!"
#endif

#include <opm/models/discretefracture/discretefracturemodel.hh>
#include <opm/models/io/dgfvanguard.hh>

#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Dnapl.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <iostream>
#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class FractureProblem;
}

BEGIN_PROPERTIES

// Create a type tag for the problem
NEW_TYPE_TAG(FractureProblem, INHERITS_FROM(DiscreteFractureModel));

// Set the grid type
SET_TYPE_PROP(
    FractureProblem, Grid,
    Dune::ALUGrid</*dim=*/2, /*dimWorld=*/2, Dune::simplex, Dune::nonconforming>);

// Set the Vanguard property
SET_TYPE_PROP(FractureProblem, Vanguard, Opm::DgfVanguard<TypeTag>);

// Set the problem property
SET_TYPE_PROP(FractureProblem, Problem, Opm::FractureProblem<TypeTag>);

// Set the wetting phase
SET_PROP(FractureProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(FractureProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::DNAPL<Scalar> > type;
};

// Set the material Law
SET_PROP(FractureProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>
    Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffectiveLaw;
    // typedef RegularizedVanGenuchten<Traits> EffectiveLaw;
    // typedef LinearMaterial<Traits> EffectiveLaw;
public:
    typedef Opm::EffToAbsLaw<EffectiveLaw> type;
};

// Enable the energy equation
SET_BOOL_PROP(FractureProblem, EnableEnergy, true);

// Set the thermal conduction law
SET_PROP(FractureProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(FractureProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Disable gravity
SET_BOOL_PROP(FractureProblem, EnableGravity, false);

// For this problem, we use constraints to specify the left boundary
SET_BOOL_PROP(FractureProblem, EnableConstraints, true);

// Set the default value for the file name of the grid
SET_STRING_PROP(FractureProblem, GridFile, "data/fracture.art.dgf");

// Set the default value for the end time
SET_SCALAR_PROP(FractureProblem, EndTime, 3e3);

// Set the default value for the initial time step size
SET_SCALAR_PROP(FractureProblem, InitialTimeStepSize, 100);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Two-phase problem which involves fractures
 *
 * The domain is initially completely saturated by the oil phase,
 * except for the left side, which is fully water saturated. Since the
 * capillary pressure in the fractures is lower than in the rock
 * matrix and the material is hydrophilic, water infiltrates through
 * the fractures and gradually pushes the oil out on the right side,
 * where the pressure is kept constant.
 */
template <class TypeTag>
class FractureProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLawParams) ThermalConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum {
        // phase indices
        wettingPhaseIdx = MaterialLaw::wettingPhaseIdx,
        nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx,

        // number of phases
        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    template <int dim>
    struct FaceLayout
    {
        bool contains(Dune::GeometryType gt)
        { return gt.dim() == dim - 1; }
    };
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, FaceLayout> FaceMapper;

    typedef Opm::FractureMapper<TypeTag> FractureMapper;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    FractureProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C

        matrixMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        matrixMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);
        fractureMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.0);
        fractureMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

#if 0 // linear
        matrixMaterialParams_.setEntryPC(0.0);
        matrixMaterialParams_.setMaxPC(2000.0);
        fractureMaterialParams_.setEntryPC(0.0);
        fractureMaterialParams_.setMaxPC(1000.0);
#endif

#if 1 // Brooks-Corey
        matrixMaterialParams_.setEntryPressure(2000);
        matrixMaterialParams_.setLambda(2.0);
        matrixMaterialParams_.setPcLowSw(1e-1);
        fractureMaterialParams_.setEntryPressure(1000);
        fractureMaterialParams_.setLambda(2.0);
        fractureMaterialParams_.setPcLowSw(5e-2);
#endif

#if 0 // van Genuchten
        matrixMaterialParams_.setVgAlpha(0.0037);
        matrixMaterialParams_.setVgN(4.7);
        fractureMaterialParams_.setVgAlpha(0.0025);
        fractureMaterialParams_.setVgN(4.7);
#endif

        matrixMaterialParams_.finalize();
        fractureMaterialParams_.finalize();

        matrixK_ = this->toDimMatrix_(1e-15);         // m^2
        fractureK_ = this->toDimMatrix_(1e5 * 1e-15); // m^2

        matrixPorosity_ = 0.10;
        fracturePorosity_ = 0.25;
        fractureWidth_ = 1e-3; // [m]

        // initialize the energy-related parameters
        initEnergyParams_(thermalConductionParams_, matrixPorosity_);
    }

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "fracture_" << Model::name();
        return oss.str();
    }

    /*!
     * \brief Called directly after the time integration.
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
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    // \}

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
    { return matrixK_; }

    /*!
     * \brief Intrinsic permeability of fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const DimMatrix& fractureIntrinsicPermeability(const Context& context OPM_UNUSED,
                                                   unsigned spaceIdx OPM_UNUSED,
                                                   unsigned timeIdx OPM_UNUSED) const
    { return fractureK_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return matrixPorosity_; }

    /*!
     * \brief The porosity inside the fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    Scalar fracturePorosity(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    { return fracturePorosity_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context OPM_UNUSED,
                                               unsigned spaceIdx OPM_UNUSED,
                                               unsigned timeIdx OPM_UNUSED) const
    { return matrixMaterialParams_; }

    /*!
     * \brief The parameters for the material law inside the fractures.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const MaterialLawParams& fractureMaterialLawParams(const Context& context OPM_UNUSED,
                                                       unsigned spaceIdx OPM_UNUSED,
                                                       unsigned timeIdx OPM_UNUSED) const
    { return fractureMaterialParams_; }

    /*!
     * \brief Returns the object representating the fracture topology.
     */
    const FractureMapper& fractureMapper() const
    { return this->simulator().vanguard().fractureMapper(); }

    /*!
     * \brief Returns the width of the fracture.
     *
     * \todo This method should get one face index instead of two
     *       vertex indices. This probably requires a new context
     *       class, though.
     *
     * \param context The execution context.
     * \param spaceIdx1 The local index of the edge's first edge.
     * \param spaceIdx2 The local index of the edge's second edge.
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar fractureWidth(const Context& context OPM_UNUSED,
                         unsigned spaceIdx1 OPM_UNUSED,
                         unsigned spaceIdx2 OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return fractureWidth_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams&
    thermalConductionLawParams(const Context& context OPM_UNUSED,
                               unsigned spaceIdx OPM_UNUSED,
                               unsigned timeIdx OPM_UNUSED) const
    { return thermalConductionParams_; }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyParams_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(pos)) {
            // on the right boundary, we impose a free-flow
            // (i.e. Dirichlet) condition
            FluidState fluidState;
            fluidState.setTemperature(temperature_);

            fluidState.setSaturation(wettingPhaseIdx, 0.0);
            fluidState.setSaturation(nonWettingPhaseIdx,
                                     1.0 - fluidState.saturation(wettingPhaseIdx));

            fluidState.setPressure(wettingPhaseIdx, 1e5);
            fluidState.setPressure(nonWettingPhaseIdx, fluidState.pressure(wettingPhaseIdx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updateAll(fluidState);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fluidState.setDensity(phaseIdx,
                                      FluidSystem::density(fluidState, paramCache, phaseIdx));
                fluidState.setViscosity(phaseIdx,
                                        FluidSystem::viscosity(fluidState, paramCache, phaseIdx));
            }

            // set a free flow (i.e. Dirichlet) boundary
            values.setFreeFlow(context, spaceIdx, timeIdx, fluidState);
        }
        else
            // for the upper, lower and left boundaries, use a no-flow
            // condition (i.e. a Neumann 0 condition)
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    // \{

    /*!
     * \copydoc FvBaseProblem::constraints
     */
    template <class Context>
    void constraints(Constraints& constraints, const Context& context,
                     unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (!onLeftBoundary_(pos))
            // only impose constraints adjacent to the left boundary
            return;

        unsigned globalIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        if (!fractureMapper().isFractureVertex(globalIdx)) {
            // do not impose constraints if the finite volume does
            // not contain fractures.
            return;
        }

        // if the current finite volume is on the left boundary
        // and features a fracture, specify the fracture fluid
        // state.
        FluidState fractureFluidState;
        fractureFluidState.setTemperature(temperature_ + 10.0);

        fractureFluidState.setSaturation(wettingPhaseIdx, 1.0);
        fractureFluidState.setSaturation(nonWettingPhaseIdx,
                                         1.0 - fractureFluidState.saturation(
                                                   wettingPhaseIdx));

        Scalar pCFracture[numPhases];
        MaterialLaw::capillaryPressures(pCFracture, fractureMaterialParams_,
                                        fractureFluidState);

        fractureFluidState.setPressure(wettingPhaseIdx, /*pressure=*/1.0e5);
        fractureFluidState.setPressure(nonWettingPhaseIdx,
                                       fractureFluidState.pressure(wettingPhaseIdx)
                                       + (pCFracture[nonWettingPhaseIdx]
                                          - pCFracture[wettingPhaseIdx]));

        constraints.setActive(true);
        constraints.assignNaiveFromFracture(fractureFluidState,
                                            matrixMaterialParams_);
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values,
                 const Context& context OPM_UNUSED,
                 unsigned spaceIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED) const
    {
        FluidState fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wettingPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(nonWettingPhaseIdx, fluidState.pressure(wettingPhaseIdx));

        fluidState.setSaturation(wettingPhaseIdx, 0.0);
        fluidState.setSaturation(nonWettingPhaseIdx,
                                 1.0 - fluidState.saturation(wettingPhaseIdx));

        values.assignNaive(fluidState);
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

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    void initEnergyParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        // assume the volumetric heat capacity of granite
        solidEnergyParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                * 2700.0); // density of granite [kg/m^3]
        solidEnergyParams_.finalize();

        Scalar lambdaGranite = 2.8; // [W / (K m)]

        // create a Fluid state which has all phases present
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(293.15);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, 1.0135e5);
        }

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updateAll(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar lambdaSaturated;
            if (FluidSystem::isLiquid(phaseIdx)) {
                Scalar lambdaFluid = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
                lambdaSaturated =
                    std::pow(lambdaGranite, (1 - poro))
                    + std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1 - poro));

            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
        }

        Scalar lambdaVac = std::pow(lambdaGranite, (1 - poro));
        params.setVacuumLambda(lambdaVac);
    }

    DimMatrix matrixK_;
    DimMatrix fractureK_;

    Scalar matrixPorosity_;
    Scalar fracturePorosity_;

    Scalar fractureWidth_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams matrixMaterialParams_;

    ThermalConductionLawParams thermalConductionParams_;
    SolidEnergyLawParams solidEnergyParams_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Opm

#endif // EWOMS_FRACTURE_PROBLEM_HH
