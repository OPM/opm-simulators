// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::RichardsLensProblem
 */
#ifndef EWOMS_RICHARDS_LENS_PROBLEM_HH
#define EWOMS_RICHARDS_LENS_PROBLEM_HH

#include <ewoms/models/richards/richardsmodel.hh>

#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidmatrixinteractions/2p/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/2p/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/2p/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/2pAdapter.hpp>

#include <dune/grid/io/file/dgfparser.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
template <class TypeTag>
class RichardsLensProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(RichardsLensProblem, VcfvRichards);

// Use 2d YaspGrid
SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::YaspGrid<2>);

// Set the physical problem to be solved
SET_PROP(RichardsLensProblem, Problem)
{ typedef Ewoms::RichardsLensProblem<TypeTag> type; };

// Set the wetting phase
SET_PROP(RichardsLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> > type;
};

// Set the material Law
SET_PROP(RichardsLensProblem, MaterialLaw)
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

// Enable gravitational acceleration
SET_BOOL_PROP(RichardsLensProblem, EnableGravity, true);

// Enable partial reassembly of the Jacobian matrix
SET_BOOL_PROP(RichardsLensProblem, EnablePartialReassemble, true);

// Enable re-use of the Jacobian matrix of the last iteration of the
// previous for the first iteration of the current time step?
SET_BOOL_PROP(RichardsLensProblem, EnableJacobianRecycling, true);

// Use forward differences to approximate the Jacobian matrix
SET_INT_PROP(RichardsLensProblem, NumericDifferenceMethod, +1);

// Set the maximum number of newton iterations of a time step
SET_INT_PROP(RichardsLensProblem, NewtonMaxIterations, 28);

// Set the "desireable" number of newton iterations of a time step
SET_INT_PROP(RichardsLensProblem, NewtonTargetIterations, 18);

// Do not write the intermediate results of the newton method
SET_BOOL_PROP(RichardsLensProblem, NewtonWriteConvergence, false);

// The default for the end time of the simulation
SET_SCALAR_PROP(RichardsLensProblem, EndTime, 3000);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(RichardsLensProblem, InitialTimeStepSize, 100);

// The default DGF file to load
SET_STRING_PROP(RichardsLensProblem, GridFile, "./grids/richardslens_24x16.dgf");
}
}

namespace Ewoms {
/*!
 * \ingroup VcfvTestProblems
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain.
 *
 * The domain is rectangular. The left and right boundaries are
 * free-flow boundaries with fixed water pressure which corrosponds to
 * a fixed saturation of \f$S_w = 0\f$ in the Richards model, the
 * bottom boundary is closed. The top boundary is also closed except
 * for an infiltration section, where water is infiltrating into an
 * initially unsaturated porous medium. This problem is very similar
 * the the \c LensProblem, with the main difference being that the domain
 * is initally fully saturated by gas instead of water and water
 * instead of a \c DNAPL infiltrates from the top.
 */
template <class TypeTag>
class RichardsLensProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // copy some indices for convenience
        pressureWIdx = Indices::pressureWIdx,
        contiWEqIdx = Indices::contiWEqIdx,

        wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex),
        nPhaseIdx = 1 - wPhaseIdx,

        numPhases = FluidSystem::numPhases,

        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };

    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    //! The parameters of the material law to be used
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    RichardsLensProblem(TimeManager &timeManager)
        : ParentType(timeManager,
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , pnRef_(1e5)
    {
        eps_ = 3e-6;
        pnRef_ = 1e5;

        lensLowerLeft_[0] = 1.0;
        lensLowerLeft_[1] = 2.0;

        lensUpperRight_[0] = 4.0;
        lensUpperRight_[1] = 3.0;

        // parameters for the Van Genuchten law
        // alpha and n
        lensMaterialParams_.setVgAlpha(0.00045);
        lensMaterialParams_.setVgN(7.3);
        outerMaterialParams_.setVgAlpha(0.0037);
        outerMaterialParams_.setVgN(4.7);

        // parameters for the linear law
        // minimum and maximum pressures
//        lensMaterialParams_.setEntryPC(0);
//        outerMaterialParams_.setEntryPC(0);
//        lensMaterialParams_.setMaxPC(0);
//        outerMaterialParams_.setMaxPC(0);

        lensK_ = this->toDimMatrix_(1e-12);
        outerK_ = this->toDimMatrix_(5e-12);
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::name
     */
    const char *name() const
    { return "lens_richards"; }

    /*!
     * \copydoc VcfvMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 10; } // -> 10°C

    /*!
     * \copydoc VcfvMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(pos))
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
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(pos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \brief Return the reference pressure [Pa] of the wetting phase.
     *
     * \copydetails Doxygen::contextParams
     */
    template <class Context>
    Scalar referencePressure(const Context &context, int spaceIdx, int timeIdx) const
    { return pnRef_; }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc VcfvProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) ||
            onRightBoundary_(pos))
        {
            const auto &materialParams = this->materialLawParams(context, spaceIdx, timeIdx);

            Scalar Sw = 0.0;
            Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            fs.setSaturation(wPhaseIdx, Sw);
            fs.setSaturation(nPhaseIdx, 1.0 - Sw);

            PhaseVector pC;
            MaterialLaw::capillaryPressures(pC, materialParams, fs);
            fs.setPressure(wPhaseIdx, pnRef_ + pC[wPhaseIdx] - pC[nPhaseIdx]);
            fs.setPressure(nPhaseIdx, pnRef_);

            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);

            // inflow of water
            massRate[contiWEqIdx] = -0.04; // kg / (m * s)

            values.setMassRate(massRate);
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
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        const auto &materialParams = this->materialLawParams(context, spaceIdx, timeIdx);

        Scalar Sw = 0.0;
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1.0 - Sw);

        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fs);
        values[pressureWIdx] = pnRef_ + (pC[wPhaseIdx] - pC[nPhaseIdx]);
    }

    /*!
     * \copydoc VcfvProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx, int timeIdx) const
    { rate = Scalar(0.0); }

    //! \}

private:
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
        Scalar lambda = (this->bboxMax()[0] - pos[0])/width;
        return onUpperBoundary_(pos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    DimMatrix lensK_;
    DimMatrix outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;

    Scalar eps_;
    Scalar pnRef_;
};
} //end namespace

#endif
