// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 *
 * \copydoc Ewoms::PowerInjectionProblem
 */
#ifndef EWOMS_POWER_INJECTION_PROBLEM_HH
#define EWOMS_POWER_INJECTION_PROBLEM_HH

#include <opm/material/fluidmatrixinteractions/2p/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/2p/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/2p/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/mp/2pAdapter.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Air.hpp>
#include <ewoms/models/immiscible/immisciblemodel.hh>
#include <ewoms/io/cubegridcreator.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <type_traits>
#include <iostream>

namespace Ewoms {

template <class TypeTag>
class PowerInjectionProblem;

//////////
// Specify the properties for the powerInjection problem
//////////
namespace Properties {

NEW_TYPE_TAG(PowerInjectionBaseProblem);

// Set the grid implementation to be used
SET_TYPE_PROP(PowerInjectionBaseProblem, Grid, Dune::YaspGrid</*dim=*/1>);

// set the GridCreator property
SET_TYPE_PROP(PowerInjectionBaseProblem, GridCreator, CubeGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(PowerInjectionBaseProblem, Problem, Ewoms::PowerInjectionProblem<TypeTag>);

// Set the wetting phase
SET_PROP(PowerInjectionBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(PowerInjectionBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::GasPhase<Scalar, Opm::Air<Scalar> > type;
};

// Set the material Law
SET_PROP(PowerInjectionBaseProblem, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::RegularizedVanGenuchten<Scalar> EffectiveLaw;
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffectiveLaw> TwoPMaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

public:
    typedef Opm::TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};

// Write out the filter velocities for this problem
SET_BOOL_PROP(PowerInjectionBaseProblem, VtkWriteFilterVelocities, true);

// Disable gravity
SET_BOOL_PROP(PowerInjectionBaseProblem, EnableGravity, false);

// define the properties specific for the power injection problem
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeX, 100.0);
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeY, 1.0);
SET_SCALAR_PROP(PowerInjectionBaseProblem, DomainSizeZ, 1.0);

SET_INT_PROP(PowerInjectionBaseProblem, CellsX, 250);
SET_INT_PROP(PowerInjectionBaseProblem, CellsY, 1);
SET_INT_PROP(PowerInjectionBaseProblem, CellsZ, 1);

// The default for the end time of the simulation
SET_SCALAR_PROP(PowerInjectionBaseProblem, EndTime, 100);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(PowerInjectionBaseProblem, InitialTimeStepSize, 1e-3);
}

/*!
 * \ingroup VcfvTestProblems
 * \brief 1D Problem with very fast injection of gas on the left.
 *
 * The velocity model is chosen in the .cc file in this problem. The
 * spatial parameters are inspired by the ones given by
 *
 * V. Jambhekar: "Forchheimer Porous-media Flow models -- Numerical
 * Investigation and Comparison with Experimental Data", Master's
 * Thesis at Institute for Modelling Hydraulic and Environmental
 * Systems, University of Stuttgart, 2011
 */
template <class TypeTag>
class PowerInjectionProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
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
    PowerInjectionProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 26.6;

        // parameters for the Van Genuchten law
        // alpha and n
        materialParams_.setVgAlpha(0.00045);
        materialParams_.setVgN(7.3);

        K_ = this->toDimMatrix_(5.73e-08); // [m^2]

        setupInitialFluidState_();
    }

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
        oss << "powerinjection_";
        if (std::is_same<typename GET_PROP_TYPE(TypeTag, VelocityModule),
                         Ewoms::VcfvDarcyVelocityModule<TypeTag> >::value)
            oss << "darcy";
        else
            oss << "forchheimer";
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
            std::cout<<"Storage: " << storage << std::endl;
        }
    }
    //! \}

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    { return K_; }

    /*!
     * \copydoc VcfvForchheimerBaseProblem::ergunCoefficient
     */
    template <class Context>
    Scalar ergunCoefficient(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.3866; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.558; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    { return materialParams_; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return temperature_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     *
     * This problem sets a very high injection rate of nitrogen on the
     * left and a free-flow boundary on the right.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiNEqIdx] = -1.00; // kg / (m^2 * s)

            // impose a forced flow boundary
            values.setMassRate(massRate);
        }
        else  {
            // free flow boundary with initial condition on the right
            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
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
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
    }


    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    void setupInitialFluidState_()
    {
        initialFluidState_.setTemperature(temperature_);

        Scalar Sw = 1.0;
        initialFluidState_.setSaturation(wPhaseIdx, Sw);
        initialFluidState_.setSaturation(nPhaseIdx, 1 - Sw);

        Scalar p = 1e5;
        initialFluidState_.setPressure(wPhaseIdx, p);
        initialFluidState_.setPressure(nPhaseIdx, p);
    }

    DimMatrix K_;
    MaterialLawParams materialParams_;

    Opm::ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;
    Scalar temperature_;
    Scalar eps_;
};

} //end namespace

#endif
