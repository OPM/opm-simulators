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
 * \copydoc Opm::LensProblem
 */
#ifndef EWOMS_LENS_PROBLEM_HH
#define EWOMS_LENS_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <opm/material/components/Dnapl.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>

#include <opm/models/common/multiphasebaseparameters.hh>
#include <opm/models/common/transfluxmodule.hh>

#include <opm/models/discretization/common/fvbaseadlocallinearizer.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

#include <opm/models/immiscible/immiscibleproperties.hh>

#include <opm/models/io/structuredgridvanguard.hh>
#include <opm/models/io/vtkmultiphasemodule.hpp>

#include <iostream>
#include <sstream>
#include <string>

namespace Opm {
template <class TypeTag>
class LensProblem;
}

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct LensBaseProblem { using InheritsFrom = std::tuple<StructuredGridVanguard>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LensBaseProblem> { using type = Opm::LensProblem<TypeTag>; };

// Use Dune-grid's YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::LensBaseProblem> { using type = Dune::YaspGrid<2>; };

// Set the wetting phase
template<class TypeTag>
struct WettingPhase<TypeTag, TTag::LensBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::SimpleH2O<Scalar> >;
};

// Set the non-wetting phase
template<class TypeTag>
struct NonwettingPhase<TypeTag, TTag::LensBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::LiquidPhase<Scalar, Opm::DNAPL<Scalar> >;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::LensBaseProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { wettingPhaseIdx = FluidSystem::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::wettingPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::nonWettingPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffectiveLaw = Opm::RegularizedVanGenuchten<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::EffToAbsLaw<EffectiveLaw>;
};

} // namespace Opm::Properties

namespace Opm::Parameters {

// define the properties specific for the lens problem
template<class Scalar>
struct LensLowerLeftX { static constexpr Scalar value = 1.0; };

template<class Scalar>
struct LensLowerLeftY { static constexpr Scalar value = 2.0;  };

template<class Scalar>
struct LensLowerLeftZ { static constexpr Scalar value = 0.0; };

template<class Scalar>
struct LensUpperRightX { static constexpr Scalar value = 4.0; };

template<class Scalar>
struct LensUpperRightY { static constexpr Scalar value = 3.0; };

template<class Scalar>
struct LensUpperRightZ { static constexpr Scalar value = 1.0; };

} // namespace Opm::Parameters

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1m, 2m) to (4m, 3m)
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
class LensProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using WettingPhase = GetPropType<TypeTag, Properties::WettingPhase>;
    using NonwettingPhase = GetPropType<TypeTag, Properties::NonwettingPhase>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wettingPhaseIdx = FluidSystem::wettingPhaseIdx,
        nonWettingPhaseIdx = FluidSystem::nonWettingPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nonWettingPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    LensProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C
        lensLowerLeft_[0] = Parameters::Get<Parameters::LensLowerLeftX<Scalar>>();
        lensLowerLeft_[1] = Parameters::Get<Parameters::LensLowerLeftY<Scalar>>();
        lensUpperRight_[0] = Parameters::Get<Parameters::LensUpperRightX<Scalar>>();
        lensUpperRight_[1] = Parameters::Get<Parameters::LensUpperRightY<Scalar>>();

        if constexpr (dim == 3) {
            lensLowerLeft_[2] = Parameters::Get<Parameters::LensLowerLeftZ<Scalar>>();
            lensUpperRight_[2] = Parameters::Get<Parameters::LensUpperRightZ<Scalar>>();
        }

        // residual saturations
        lensMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.18);
        lensMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);
        outerMaterialParams_.setResidualSaturation(wettingPhaseIdx, 0.05);
        outerMaterialParams_.setResidualSaturation(nonWettingPhaseIdx, 0.0);

        // parameters for the Van Genuchten law: alpha and n
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
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        Parameters::Register<Parameters::LensLowerLeftX<Scalar>>
            ("The x-coordinate of the lens' lower-left corner [m].");
        Parameters::Register<Parameters::LensLowerLeftY<Scalar>>
            ("The y-coordinate of the lens' lower-left corner [m].");
        Parameters::Register<Parameters::LensUpperRightX<Scalar>>
            ("The x-coordinate of the lens' upper-right corner [m].");
        Parameters::Register<Parameters::LensUpperRightY<Scalar>>
            ("The y-coordinate of the lens' upper-right corner [m].");

        if constexpr (dim == 3) {
            Parameters::Register<Parameters::LensLowerLeftZ<Scalar>>
                ("The z-coordinate of the lens' lower-left corner [m].");
            Parameters::Register<Parameters::LensUpperRightZ<Scalar>>
                ("The z-coordinate of the lens' upper-right corner [m].");
        }

        Parameters::SetDefault<Parameters::CellsX>(48);
        Parameters::SetDefault<Parameters::CellsY>(32);
        Parameters::SetDefault<Parameters::DomainSizeX<Scalar>>(6.0);
        Parameters::SetDefault<Parameters::DomainSizeY<Scalar>>(4.0);

        if constexpr (dim == 3) {
            Parameters::SetDefault<Parameters::CellsZ>(16);
            Parameters::SetDefault<Parameters::DomainSizeZ<Scalar>>(1.0);
        }

        // Use forward differences
        using LLS = GetPropType<TypeTag, Properties::LocalLinearizerSplice>;
        constexpr bool useFD = std::is_same_v<LLS, Properties::TTag::FiniteDifferenceLocalLinearizer>;
        if constexpr (useFD) {
            Parameters::SetDefault<Parameters::NumericDifferenceMethod>(+1);
        }

        Parameters::SetDefault<Parameters::EndTime<Scalar>>(30e3);
        Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);
        Parameters::SetDefault<Parameters::EnableStorageCache>(true);
        Parameters::SetDefault<Parameters::InitialTimeStepSize<Scalar>>(250.0);
        Parameters::SetDefault<Parameters::VtkWriteIntrinsicPermeabilities>(true);
        Parameters::SetDefault<Parameters::EnableGravity>(true);
    }

    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription()
    {
        std::string thermal = "isothermal";
        constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
        if constexpr (enableEnergy)
            thermal = "non-isothermal";

        std::string deriv = "finite difference";
        using LLS = GetPropType<TypeTag, Properties::LocalLinearizerSplice>;
        constexpr bool useAutoDiff = std::is_same_v<LLS, Properties::TTag::AutoDiffLocalLinearizer>;
        if constexpr (useAutoDiff) {
            deriv = "automatic differentiation";
        }

        std::string disc = "vertex centered finite volume";
        using D = GetPropType<TypeTag, Properties::Discretization>;
        constexpr bool useEcfv = std::is_same<D, Opm::EcfvDiscretization<TypeTag>>::value;
        if constexpr (useEcfv)
            disc = "element centered finite volume";

        return std::string("")+
            "Ground remediation problem where a dense oil infiltrates "+
            "an aquifer with an embedded low-permability lens. " +
            "This is the binary for the "+thermal+" variant using "+deriv+
            "and the "+disc+" discretization";
    }

    /*!
     * \name Soil parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);

        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& /*context*/,
                    unsigned /*spaceIdx*/,
                    unsigned /*timeIdx*/) const
    { return 0.4; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& globalPos = context.pos(spaceIdx, timeIdx);

        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& /*context*/,
                       unsigned /*spaceIdx*/,
                       unsigned /*timeIdx*/) const
    { return temperature_; }

    //! \}

    /*!
     * \name Auxiliary methods
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        using LLS = GetPropType<TypeTag, Properties::LocalLinearizerSplice>;

        constexpr bool useAutoDiff = std::is_same_v<LLS, Properties::TTag::AutoDiffLocalLinearizer>;

        using FM = GetPropType<TypeTag, Properties::FluxModule>;
        constexpr bool useTrans = std::is_same_v<FM, Opm::TransFluxModule<TypeTag>>;

        std::ostringstream oss;
        oss << "lens_" << Model::name()
            << "_" << Model::discretizationName()
            << "_" << (useAutoDiff?"ad":"fd");
        if (useTrans)
            oss << "_trans";

        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::beginTimeStep
     */
    void beginTimeStep()
    { }

    /*!
     * \copydoc FvBaseProblem::beginIteration
     */
    void beginIteration()
    { }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
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

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            // free flow boundary. we assume incompressible fluids
            Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/Scalar(1e5));
            Scalar densityN = NonwettingPhase::density(temperature_, /*pressure=*/Scalar(1e5));

            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;

            // set wetting phase pressure and saturation
            if (onLeftBoundary_(pos)) {
                Scalar height = this->boundingBoxMax()[1] - this->boundingBoxMin()[1];
                Scalar depth = this->boundingBoxMax()[1] - pos[1];
                Scalar alpha = (1 + 1.5 / height);

                // hydrostatic pressure scaled by alpha
                pw = 1e5 - alpha * densityW * this->gravity()[1] * depth;
                Sw = 1.0;
            }
            else {
                Scalar depth = this->boundingBoxMax()[1] - pos[1];

                // hydrostatic pressure
                pw = 1e5 - densityW * this->gravity()[1] * depth;
                Sw = 1.0;
            }

            // specify a full fluid state using pw and Sw
            const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);

            Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                                      /*storeEnthalpy=*/false> fs;
            fs.setSaturation(wettingPhaseIdx, Sw);
            fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);
            fs.setTemperature(T);

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wettingPhaseIdx, pw);
            fs.setPressure(nonWettingPhaseIdx, pw + pC[nonWettingPhaseIdx] - pC[wettingPhaseIdx]);

            fs.setDensity(wettingPhaseIdx, densityW);
            fs.setDensity(nonWettingPhaseIdx, densityN);

            fs.setViscosity(wettingPhaseIdx, WettingPhase::viscosity(temperature_, fs.pressure(wettingPhaseIdx)));
            fs.setViscosity(nonWettingPhaseIdx, NonwettingPhase::viscosity(temperature_, fs.pressure(nonWettingPhaseIdx)));

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
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        Scalar depth = this->boundingBoxMax()[1] - pos[1];

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wettingPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = 1.0;
        fs.setSaturation(wettingPhaseIdx, Sw);
        fs.setSaturation(nonWettingPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fs, wettingPhaseIdx);
        Scalar densityW = FluidSystem::density(fs, paramCache, wettingPhaseIdx);

        // hydrostatic pressure (assuming incompressibility)
        Scalar pw = 1e5 - densityW * this->gravity()[1] * depth;

        // calculate the capillary pressure
        const MaterialLawParams& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // make a full fluid state
        fs.setPressure(wettingPhaseIdx, pw);
        fs.setPressure(nonWettingPhaseIdx, pw + (pC[wettingPhaseIdx] - pC[nonWettingPhaseIdx]));

        // assign the primary variables
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
                const Context& /*context*/,
                unsigned /*spaceIdx*/,
                unsigned /*timeIdx*/) const
    { rate = Scalar(0.0); }

    //! \}

private:
    bool isInLens_(const GlobalPosition& pos) const
    {
        for (unsigned i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] - eps_ || pos[i] > lensUpperRight_[i]
                                                              + eps_)
                return false;
        }
        return true;
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < this->boundingBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[1] < this->boundingBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[1] > this->boundingBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    {
        Scalar width = this->boundingBoxMax()[0] - this->boundingBoxMin()[0];
        Scalar lambda = (this->boundingBoxMax()[0] - pos[0]) / width;
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

} // namespace Opm

#endif
