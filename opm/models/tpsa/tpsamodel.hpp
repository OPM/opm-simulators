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
#ifndef TPSA_MODEL_HPP
#define TPSA_MODEL_HPP

#include <dune/grid/common/gridenums.hh>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <array>
#include <memory>


namespace Opm {

/*!
* \brief TPSA geomechanics model
*
* Solves the (linear) elasticity equations using Newton method.
*/
template <class TypeTag>
class TpsaModel
{
    using DofMapper = GetPropType<TypeTag, Properties::DofMapper>;
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVectorTPSA>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Indices = GetPropType<TypeTag, Properties::IndicesTPSA>;
    using Linearizer = GetPropType<TypeTag, Properties::LinearizerTPSA>;
    using NewtonMethod = GetPropType<TypeTag, Properties::NewtonMethodTPSA>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariablesTPSA>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVectorTPSA>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    enum { dimWorld = GridView::dimensionworld };
    enum { historySize = getPropValue<TypeTag, Properties::SolutionHistorySizeTPSA>() };
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };

    enum { disp0Idx = Indices::disp0Idx };
    enum { rot0Idx = Indices::rot0Idx };
    enum { solidPres0Idx = Indices::solidPres0Idx };

    using MaterialState = MaterialStateTPSA<Evaluation>;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using SymTensor = Dune::FieldVector<Scalar, 6>;

public:
    /*!
    * \brief Small block vector wrapper class for model solutions
    */
    class TpsaBlockVectorWrapper
    {
    protected:
        SolutionVector blockVector_;
    public:
        /*!
        * \brief Constructor
        *
        * \param
        * \param size Size of block vector
        */
        TpsaBlockVectorWrapper(const std::string&, const std::size_t size)
            : blockVector_(size)
        {}

        /*!
        * \brief Default constructor
        */
        TpsaBlockVectorWrapper() = default;

        /*!
        * \brief Test function for serialization
        */
        static TpsaBlockVectorWrapper serializationTestObject()
        {
            TpsaBlockVectorWrapper result("dummy", 3);
            result.blockVector_[0] = 1.0;
            result.blockVector_[1] = 2.0;
            result.blockVector_[2] = 3.0;

            return result;
        }

        /*!
        * \brief Get reference of block vector
        *
        * \returns Block vector
        */
        SolutionVector& blockVector()
        { return blockVector_; }

        /*!
        * \brief Get const reference of block vector
        *
        * \returns Block vector
        */
        const SolutionVector& blockVector() const
        { return blockVector_; }

        /*!
        * \brief Check if incoming block vector is the same as current
        *
        * \param wrapper Block vector to check
        * \returns Boolean indicating if wrapper is equal to current block vector
        */
        bool operator==(const TpsaBlockVectorWrapper& wrapper) const
        {
            return std::equal(this->blockVector_.begin(), this->blockVector_.end(),
                              wrapper.blockVector_.begin(), wrapper.blockVector_.end());
        }

        /*!
        * \brief Serializing operation
        *
        * \param serializer Reference to serializer object
        */
        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(blockVector_);
        }
    };

    // ///
    // Public functions
    // ///
    /*!
    * \brief Constructor
    *
    * \param simulator Simulator object
    */
    explicit TpsaModel(Simulator& simulator)
        : linearizer_(std::make_unique<Linearizer>())
        , newtonMethod_(simulator)
        , simulator_(simulator)
        , element_chunks_(simulator.gridView(), Dune::Partitions::all, ThreadManager::maxThreads())
    {
        // Initialize equation weights to 1.0
        eqWeights_.resize(numEq, 1.0);

        // Initialize historic solution vectors
        // OBS: need at least history size = 2, due to time-derivative of solid-pressure in Flow coupling term
        const std::size_t numDof = simulator_.model().numGridDof();
        for (unsigned timeIdx = 0; timeIdx < historySize; ++timeIdx) {
            solution_[timeIdx] = std::make_unique<TpsaBlockVectorWrapper>("solution", numDof);
        }
    }

    /*!
    * \brief Initialize TPSA model
    */
    void finishInit()
    {
        // Initialize the linearizer
        linearizer_->init(simulator_);

        // Resize material state vector
        resizeMaterialState_();
    }

    /*!
    * \brief Register runtime parameters
    */
    static void registerParameters()
    {
        // Newton method parameters
        NewtonMethod::registerParameters();
    }

    /*!
    * \brief Prepare TPSA model for coupled Flow-TPSA scheme
    */
    void prepareTPSA()
    {
        // Update historic solution
        solution(/*timeIdx=*/1) = solution(/*timeIdx=*/0);
    }

    /*!
    * \brief Sync primary variables in overlapping cells
    *
    * Copied code from EcfvDiscretization::syncOverlap() to sync TPSA primary variables
    */
    void syncOverlap()
    {
        // Syncronize the solution on the ghost and overlap elements
        using GhostSyncHandle = GridCommHandleGhostSync<PrimaryVariables,
                                                        SolutionVector,
                                                        DofMapper,
                                                        /*commCodim=*/0>;

        auto ghostSync = GhostSyncHandle(solution(/*timeIdx=*/0),
                                         simulator_.model().dofMapper());

        simulator_.gridView().communicate(ghostSync,
                                          Dune::InteriorBorder_All_Interface,
                                          Dune::ForwardCommunication);
    }

    /*!
    * \brief Update material state for all cells
    *
    * \param timeIdx Time index
    *
    * \note Cached material state not implemented yet!
    */
    void updateMaterialState(const unsigned /*timeIdx*/)
    {
        // Loop over all elements chuncks and update material state from current solution
        const auto& elementMapper = simulator_.model().elementMapper();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (const auto& chunk : element_chunks_) {
            for (const auto& elem : chunk) {
                const unsigned globalIdx = elementMapper.index(elem);
                auto& currSol = solution(/*timeIdx=*/0)[globalIdx];
                setMaterialState_(globalIdx, /*timeIdx=*/0, currSol);
            }
        }
    }

    // ///
    // Public get and set functions
    // ///
    /*!
    * \brief Return the linearizer
    *
    * \returns Reference to linearizer
    */
    const Linearizer& linearizer() const
    {
        return *linearizer_;
    }

    /*!
    * \brief Return the linearizer
    *
    * \returns Reference to linearizer
    */
    Linearizer& linearizer()
    {
        return *linearizer_;
    }

    /*!
    * \brief Return the Newton method
    *
    * \returns Reference to Newton method object
    */
    const NewtonMethod& newtonMethod() const
    {
        return newtonMethod_;
    }

    /*!
    * \brief Return the Newton method
    *
    * \returns Reference to Newton method object
    */
    NewtonMethod& newtonMethod()
    {
        return newtonMethod_;
    }

    /*!
    * \brief Get reference to history solution vector
    *
    * \param timeIdx Time index
    * \returns Reference to solution vector at time step timeIdx
    */
    const SolutionVector& solution(unsigned timeIdx) const
    {
        return solution_[timeIdx]->blockVector();
    }

    /*!
    * \brief Get reference to history solution vector
    *
    * \param timeIdx Time index
    * \returns Reference to solution vector at time step timeIdx
    */
    SolutionVector& solution(unsigned timeIdx)
    {
        return solution_[timeIdx]->blockVector();
    }

    /*!
    * \brief Return number of degrees of freedom in the grid from the Flow model
    * \returns Number of grid DOFs
    */
    std::size_t numGridDof() const
    {
        return simulator_.model().numGridDof();
    }

    /*!
    * \brief Return the total number of degrees of freedom
    * \returns Number of grid DOFs + auxillary DOFs
    */
    std::size_t numTotalDof() const
    {
        return numGridDof() + numAuxiliaryDof();
    }

    /*!
    * \brief Return the total grid volume from the Flow model
    *
    * \param globalIdx Cell index
    * \returns Grid volume of grid cell
    */
    Scalar dofTotalVolume(unsigned globalIdx) const
    {
        return simulator_.model().dofTotalVolume(globalIdx);
    }

    /*!
    * \brief Return equation weights
    *
    * \param dofIdx Degree-of-freedom index
    * \param eqIdx Equation index
    * \returns Weight for equation
    */
    Scalar eqWeight(unsigned /*dofIdx*/, unsigned eqIdx) const
    {
        return eqWeights_[eqIdx];
    }

    /*!
    * \brief Set weights for equation
    *
    * \param eqIdx Equation index
    * \param value Weight value to set
    */
    void setEqWeight(unsigned eqIdx, Scalar value)
    {
        eqWeights_[eqIdx] = value;
    }

    /*!
    * \brief Return number of auxillary modules
    * \returns No. of auxillary modules
    */
    std::size_t numAuxiliaryModules() const
    {
        return 0;
    }

    /*!
    * \brief Return number of auxillary degrees of freedom
    * \returns No. of auxillary DOFs
    */
    std::size_t numAuxiliaryDof() const
    {
        return 0;
    }

    /*!
    * \brief Return current material state
    *
    * \param globalIdx Cell index
    * \param timeIdx Time index
    * \returns Material state container for grid cell
    *
    * \note Cached material state not implemented yet!
    */
    const MaterialState& materialState(const unsigned globalIdx, unsigned /*timeIdx*/) const
    {
        return materialState_[globalIdx];
    }

    /*!
    * \brief Output displacement vector
    *
    * \param globalIdx Cell index
    * \param with_fracture Boolean to activate fracture output
    * \returns Displacement vector at grid cell
    *
    * \note Needed in OutputBlackOilModule!
    */
    DimVector disp(const unsigned globalIdx, const bool /*with_fracture*/) const
    {
        DimVector d;
        for (std::size_t i = 0; i < 3; ++i) {
            d[i] = decay<Scalar>(materialState_[globalIdx].displacement(i));
        }
        return d;
    }

    /*!
    * \brief Output (del?)stress tensor
    *
    * \param globalIdx Cell index
    * \returns (Del?)stress tensor (Voigt notation) at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    SymTensor delstress(const unsigned /*globalIdx*/) const
    {
        SymTensor val;
        return val;
    }

    /*!
    * \brief Output fracture stress tensor
    *
    * \param globalIdx Cell index
    * \returns Fracture stress tensor (Voigt notation) at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    SymTensor fractureStress(const unsigned /*globalIdx*/) const
    {
        SymTensor val;
        return val;
    }

    /*!
    * \brief Output linear stress tensor
    *
    * \param globalIdx Cell index
    * \returns Linear stress tensor (Voigt notation) at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    SymTensor linstress(const unsigned /*globalIdx*/) const
    {
        SymTensor val;
        return val;
    }

    /*!
    * \brief Output stress tensor
    *
    * \param globalIdx Cell index
    * \param with_fracture Boolean to activate fracture output
    * \returns Stress tensor (Voigt notation) at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    SymTensor stress(const unsigned /*globalIdx*/, const bool /*with_fracture*/) const
    {
        SymTensor val;
        return val;
    }

    /*!
    * \brief Output strain tensor
    *
    * \param globalIdx Cell index
    * \param with_fracture Boolean to activate fracture output
    * \returns Strain tensor (Voigt notation) at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    SymTensor strain(const unsigned /*globalIdx*/, const bool /*with_fracture*/) const
    {
        SymTensor val;
        return val;
    }

    /*!
    * \brief Output potential forces
    *
    * \param globalIdx Cell index
    * \returns Potential forces at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    Scalar mechPotentialForce(unsigned /*globalIdx*/) const
    {
        return Scalar(0.0);
    }

    /*!
    * \brief Output potential pressure forces
    *
    * \param globalIdx Cell index
    * \returns Potential pressure forces at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    Scalar mechPotentialPressForce(unsigned /*globalIdx*/) const
    {
        return Scalar(0.0);
    }

    /*!
    * \brief Output potential temparature forces
    *
    * \param globalIdx Cell index
    * \returns Potential temperature forces at grid cell
    *
    * \note Needed in OutputBlackOilModule, but zero for now!
    */
    Scalar mechPotentialTempForce(unsigned /*globalIdx*/) const
    {
        return Scalar(0.0);
    }

protected:
    // ///
    // Protected functions
    // ///
    /*!
    * \brief Resize material state vector
    *
    * \note Cached material state not implemented yet!
    */
    void resizeMaterialState_()
    {
        const std::size_t numDof = simulator_.model().numGridDof();
        materialState_.resize(numDof);
    }

private:
    // ///
    // Private functions
    // ///
    /*!
    * \brief Set material from another
    *
    * \param globalIdx Cell index
    * \param timeIdx Time index
    * \param values Primary variable to insert into material state
    *
    * \note Cached material state not implemented yet!
    */
    void setMaterialState_(const unsigned globalIdx, const unsigned /*timeIdx*/, PrimaryVariables& values)
    {
        auto& dofMaterialState = materialState_[globalIdx];
            for (unsigned dirIdx = 0; dirIdx < 3; ++dirIdx) {
                dofMaterialState.setDisplacement(dirIdx, values.makeEvaluation(disp0Idx + dirIdx, 0));
                dofMaterialState.setRotation(dirIdx, values.makeEvaluation(rot0Idx + dirIdx, 0));
            }
            dofMaterialState.setSolidPressure(values.makeEvaluation(solidPres0Idx, 0));
    }

    std::unique_ptr<Linearizer> linearizer_;
    NewtonMethod newtonMethod_;
    Simulator& simulator_;
    ElementChunks<GridView, Dune::Partitions::All> element_chunks_;

    std::array<std::unique_ptr<TpsaBlockVectorWrapper>, historySize> solution_;
    std::vector<Scalar> eqWeights_;
    std::vector<MaterialState> materialState_;
};  // class TpsaModel

}  // namespace Opm


#endif
