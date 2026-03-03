// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef FLOW_PROBLEM_TPSA_HPP
#define FLOW_PROBLEM_TPSA_HPP

#include <dune/common/fvector.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/materialstates/MaterialStateTPSA.hpp>

#include <opm/models/io/vtktpsamodule.hpp>
#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/tpsa/tpsamodel.hpp>
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/FacePropertiesTPSA.hpp>
#include <opm/simulators/flow/FlowProblemBlackoil.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <fmt/format.h>


namespace Opm {

/*!
* \brief Problem for Flow-TPSA coupled simulations
*/
template <class TypeTag>
class FlowProblemTPSA : public FlowProblemBlackoil<TypeTag>
{
public:
    using ParentType = FlowProblemBlackoil<TypeTag>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GeomechModel = GetPropType<TypeTag, Properties::ModelTPSA>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::IndicesTPSA>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum { dimWorld = GridView::dimensionworld };
    enum { enableMech = getPropValue<TypeTag, Properties::EnableMech>() };
    enum { historySize = getPropValue<TypeTag, Properties::SolutionHistorySizeTPSA>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };
    enum { numPhases = FluidSystem::numPhases };

    enum { contiRotEqIdx = Indices::contiRotEqIdx };
    enum { contiSolidPresEqIdx = Indices::contiSolidPresEqIdx };
    enum { solidPres0Idx = Indices::solidPres0Idx };

    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using FaceProperties = FacePropertiesTPSA<Grid, GridView, ElementMapper, CartesianIndexMapper, Scalar>;
    using InitialMaterialState = MaterialStateTPSA<Scalar>;
    using Toolbox = MathToolbox<Evaluation>;

    // ///
    // Public functions
    // ///
    /*!
    * \brief Constructor
    *
    * \param simulator Reference to simulator object
    */
    FlowProblemTPSA(Simulator& simulator)
        : ParentType(simulator)
        , faceProps_(simulator.vanguard().eclState(),
                     simulator.vanguard().gridView(),
                     simulator.vanguard().cartesianIndexMapper(),
                     simulator.vanguard().grid(),
                     simulator.vanguard().cellCentroids())
        , geoMechModel_(simulator)
    {
        if constexpr(enableMech) {
            // Add VTK TPSA to output module
            this->model().addOutputModule(std::make_unique<VtkTpsaModule<TypeTag>>(simulator));

            // Sanity check
            const auto& mechSolver = simulator.vanguard().eclState().runspec().mechSolver();
            if (!mechSolver.tpsa()) {
                std::string msg = "Simulator with Tpsa-geomechanics enabled compile time, but deck does not contain "
                                  "TPSA keyword!";
                OpmLog::error(msg);
                throw std::runtime_error(msg);
            }
        }
        else {
            // Sanity check
            const auto& mechSolver = simulator.vanguard().eclState().runspec().mechSolver();
            if (mechSolver.tpsa()) {
                std::string msg = "TPSA keyword in deck, but Tpsa-geomechanics disabled compile-time!";
                OpmLog::error(msg);
                throw std::runtime_error(msg);
            }
        }
    }

    /*!
    * \brief Register runtime parameters
    */
    static void registerParameters()
    {
        // Register parameters for parent class
        ParentType::registerParameters();

        // Geomech model parameters
        GeomechModel::registerParameters();

        // VTK output parameters
        VtkTpsaModule<TypeTag>::registerParameters();
    }

    /*!
    * \brief Initialize the problem
    */
    void finishInit()
    {
        // FlowProblemBlackoil::finishInit()
        ParentType::finishInit();

        // Read initial conditions and set material state
        readInitalConditionsTPSA_();

        // Calculate face properties
        faceProps_.finishInit();

        // Set equation weights
        computeAndSetEqWeights_();

        // Initialize the TPSA model
        geoMechModel_.finishInit();
    }

    /*!
    * \brief Set initial solution for the problem
    *
    * This function is a combination FvBaseDiscretization::applyInitialSolution and FlowProblemBlackoil::initial()
    */
    void initialSolutionApplied()
    {
        // Set up initial solution for the Flow model
        ParentType::initialSolutionApplied();

        // Initialize soultions as zero
        auto& uCur = geoMechModel_.solution(/*timeIdx=*/0);
        uCur = Scalar(0.0);

        // Loop through grid and set initial solution from material state
        ElementContext elemCtx(this->simulator());
        for (const auto& elem : elements(this->gridView())) {
            // Ignore everything which is not in the interior if the current process' piece of the grid
            if (elem.partitionType() != Dune::InteriorEntity) {
                continue;
            }

            // Loop over sub control volumes and set initial solutions
            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                auto& priVars = uCur[globalIdx];
                priVars.assignNaive(initialMaterialState_[globalIdx]);
                priVars.checkDefined();
            }
        }

        // synchronize the ghost/overlapping DOFs (if necessary)
        geoMechModel_.syncOverlap();

        // Set history solutions to the initial solution.
        for (unsigned timeIdx = 1; timeIdx < historySize; ++timeIdx) {
            geoMechModel_.solution(timeIdx) = geoMechModel_.solution(/*timeIdx=*/0);
        }

        // Set material state
        geoMechModel_.updateMaterialState(/*timeIdx=*/0);
    }

    /*!
    * \brief Compute weights to rescale the TPSA equations
    */
    void computeAndSetEqWeights_()
    {
        // Average shear modulus over domain
        Scalar avgSmodulus = 0.0;
        const auto& gridView = this->gridView();
        ElementContext elemCtx(this->simulator());
        for(const auto& elem: elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            avgSmodulus += this->shearModulus(elemIdx);
        }
        std::size_t numDof = this->model().numGridDof();
        const auto& comm = this->simulator().vanguard().grid().comm();
        comm.sum(avgSmodulus);
        Scalar numTotalDof = comm.sum(numDof);
        avgSmodulus /= numTotalDof;
        avgSmodulus = std::sqrt(avgSmodulus);

        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (eqIdx < contiRotEqIdx) {
                geoMechModel_.setEqWeight(eqIdx, 1 / avgSmodulus);
            }
            else {
                geoMechModel_.setEqWeight(eqIdx, avgSmodulus);
            }
        }
    }

    /*!
    * \brief Called by the simulator before each time integration.
    */
    void beginTimeStep() override
    {
        // Call parent class beginTimeStep()
        ParentType::beginTimeStep();

        // Update mechanics boundary conditions.
        // NOTE: Flow boundary conditions should be updated in ParentType::beginTimeStep()
        if (this->nonTrivialBoundaryConditions()) {
            geoMechModel_.linearizer().updateBoundaryConditionData();
        }
    }

    /*!
    * \brief Organize mechanics boundary conditions
    *
    * \param globalSpaceIdx Cell index
    * \param directionId Direction id
    * \returns A pair of BCMECHType and displacement vector
    *
    * Output from this function is used in LocalResidual::computeBoundaryTerm
    *
    * \note Only BCMECHTYPE = FREE and NONE implemented. FIXED will/should throw an error when computed in local
    * residual!
    */
    std::pair<BCMECHType, Dune::FieldVector<Evaluation, 3>>
    mechBoundaryCondition(const unsigned int globalSpaceIdx, const int directionId)
    {
        // Default boundary conditions if BCCON/BCPROP not defined
        if (!this->nonTrivialBoundaryConditions()) {
            return { BCMECHType::NONE, Dune::FieldVector<Evaluation, 3>{0.0, 0.0, 0.0} };
        }

        // Default for BCPROP index = 0 or no BCPROP defined at current episode
        FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);
        const auto& schedule = this->simulator().vanguard().schedule();
        if (this->bcindex_(dir)[globalSpaceIdx] == 0 || schedule[this->episodeIndex()].bcprop.size() == 0) {
            return { BCMECHType::NONE, Dune::FieldVector<Evaluation, 3>{0.0, 0.0, 0.0} };
        }

        // Get current BC
        const auto& bc = schedule[this->episodeIndex()].bcprop[this->bcindex_(dir)[globalSpaceIdx]];
        if (bc.bcmechtype == BCMECHType::FREE) {
            return { BCMECHType::FREE, Dune::FieldVector<Evaluation, 3>{0.0, 0.0, 0.0} };
        }
        else {
            return { bc.bcmechtype, Dune::FieldVector<Evaluation, 3>{0.0, 0.0, 0.0} };
        }
    }

    /*!
    * \brief Set mechanics source term, in particular coupling terms
    *
    * \param sourceTerm Source term vector
    * \param globalSpaceIdx Cell index
    * \param timeIdx Time index
    *
    * This term requires that intensive quantities are updated!
    */
    void tpsaSource(Dune::FieldVector<Evaluation, numEq>& sourceTerm,
                    unsigned globalSpaceIdx,
                    unsigned timeIdx)
    {
        sourceTerm = 0.0;

        // Coupling term Flow -> TPSA
        const auto biot = this->biotCoeff(globalSpaceIdx);
        const auto lameParam = this->lame(globalSpaceIdx);

        const auto& iq = this->model().intensiveQuantities(globalSpaceIdx, timeIdx);
        const auto& fs = iq.fluidState();
        const auto pres = decay<Scalar>(fs.pressure(this->refPressurePhaseIdx_()));
        const auto initPres = this->initialFluidState(globalSpaceIdx).pressure(this->refPressurePhaseIdx_());

        auto sourceFromFlow = -biot / lameParam * (pres - initPres);
        sourceTerm[contiSolidPresEqIdx] += sourceFromFlow;
    }

    /*!
    * \brief Pore volume change due to geomechanics
    *
    * \param elementIdx Cell index
    * \param timeIdx Time index
    * \returns Pore volume change (dimensionless)
    *
    * \note This is the coupling term to Flow
    */
    Scalar rockMechPoroChange(unsigned elementIdx, unsigned timeIdx) const
    {
        // TODO: get timeIdx=1 solid pressure from a cached materialState (or intensiveQuantities) if/when implemented
        assert (timeIdx <= historySize);
        const auto solidPres = (timeIdx == 0) ?
           decay<Scalar>( geoMechModel_.materialState(elementIdx, /*timeIdx=*/timeIdx).solidPressure()) :
           geoMechModel_.solution(/*timeIdx=*/timeIdx)[elementIdx][solidPres0Idx];
        const auto biot = this->biotCoeff(elementIdx);
        const auto lameParam = this->lame(elementIdx);

        return biot / lameParam * solidPres;
    }

    // ///
    // Public get functions
    // ///
    /*!
    * \brief Direct access to average (half-)weight at interface between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    * \returns Weight average
    */
    Scalar weightAverage(unsigned globalElemIdxIn, unsigned globalElemIdxOut)
    {
        return faceProps_.weightAverage(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    * \returns Weight average at boundary
    */
    Scalar weightAverageBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.weightAverageBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to product of weights at interface between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    * \returns Weight product
    */
    Scalar weightProduct(unsigned globalElemIdxIn, unsigned globalElemIdxOut) const
    {
        return faceProps_.weightProduct(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    * \returns Normal distance
    */
    Scalar normalDistance(unsigned globalElemIdxIn, unsigned globalElemIdxOut) const
    {
        return faceProps_.normalDistance(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to normal distance at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    * \returns Normal distance to boundary
    */
    Scalar normalDistanceBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.normalDistanceBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to face normal between two elements
    *
    * \param globalElemIdxIn Inside cell index
    * \param globalElemIdxOut Outside cell index
    * \returns Face normal
    */
    DimVector cellFaceNormal(unsigned globalElemIdxIn, unsigned globalElemIdxOut)
    {
        return faceProps_.cellFaceNormal(globalElemIdxIn, globalElemIdxOut);
    }

    /*!
    * \brief Direct access to face normal at the boundary
    *
    * \param globalElemIdxIn Inside cell index
    * \param boundaryFaceIdx Boundary (local) face index
    * \returns Face normal on boundary
    */
    const DimVector& cellFaceNormalBoundary(unsigned globalElemIdxIn, unsigned boundaryFaceIdx) const
    {
        return faceProps_.cellFaceNormalBoundary(globalElemIdxIn, boundaryFaceIdx);
    }

    /*!
    * \brief Direct access to shear modulus in an element
    *
    * \param globalElemIdx Cell index
    * \returns Shear modulus
    */
    Scalar shearModulus(unsigned globalElemIdx) const
    {
        return faceProps_.shearModulus(globalElemIdx);
    }

    /*!
    * \brief Flow-TPSA lagged coupling scheme activated?
    *
    * \returns Bool indicating lagged coupling
    */
    bool laggedScheme() const
    {
        const auto& mechSolver = this->simulator().vanguard().eclState().runspec().mechSolver();
        return mechSolver.laggedScheme();
    }

    /*!
    * \brief Flow-TPSA fixed-stress coupling scheme activated?
    *
    * \returns Bool indicating fixed-stress coupling
    */
    bool fixedStressScheme() const
    {
        const auto& mechSolver = this->simulator().vanguard().eclState().runspec().mechSolver();
        return mechSolver.fixedStressScheme();
    }

    /*!
    * \brief Get TPSA model
    *
    * \returns Reference to geomechanics model
    */
    const GeomechModel& geoMechModel() const
    {
        return geoMechModel_;
    }

    /*!
    * \brief Get TPSA model
    *
    * \returns Reference to geomechanics model
    */
    GeomechModel& geoMechModel()
    {
        return geoMechModel_;
    }

    /*!
    * \brief Get fixed-stress iteration parameters
    *
    * \returns Pair with min/max fixed-stress iterations
    */
    std::pair<int, int> fixedStressParameters() const
    {
        const auto& mechSolver = this->simulator().vanguard().eclState().runspec().mechSolver();
        return std::make_pair(mechSolver.fixedStressMinIter(), mechSolver.fixedStressMaxIter());
    }

protected:
    // ///
    // Protected functions
    // ///
    /*!
    * \brief Read initial conditions and generate material state for TPSA model
    */
    void readInitalConditionsTPSA_()
    {
        // ///
        // OBS: No equilibration keywords (e.g., STREQUIL) implemented yet!
        // ///

        // Set all initial material state variables to zero
        std::size_t numDof = this->model().numGridDof();
        initialMaterialState_.resize(numDof);
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofMaterialState = initialMaterialState_[dofIdx];
            for (unsigned dirIdx = 0; dirIdx < 3; ++dirIdx) {
                dofMaterialState.setDisplacement(dirIdx, 0.0);
                dofMaterialState.setRotation(dirIdx, 0.0);
            }
            dofMaterialState.setSolidPressure(0.0);
        }
    }

private:
    FaceProperties faceProps_;
    GeomechModel geoMechModel_;

    std::vector<Scalar> biotcoeff_;
    std::vector<InitialMaterialState> initialMaterialState_;
};  // class FlowProblemTPSA

}  // namespace Opm

#endif
