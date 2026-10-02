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
 * \copydoc Opm::BlackOilNewtonMethod
 */
#ifndef OPM_BLACK_OIL_NEWTON_METHOD_HPP
#define OPM_BLACK_OIL_NEWTON_METHOD_HPP

#include <opm/common/Exceptions.hpp>

#include <opm/models/blackoil/blackoilmodules.hpp>
#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>
#include <opm/models/blackoil/blackoilnewtonupdate.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/nonlinear/newtonmethod.hh>

#include <opm/models/utils/signum.hh>

#include <opm/material/common/Valgrind.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    static constexpr bool enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>();
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag, enableBioeffects>;

    static const unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

public:
    explicit BlackOilNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        bparams_.read();
    }

    /*!
     * \copydoc NewtonMethod::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        wasSwitched_.resize(this->model().numTotalDof(), false);
    }

    void resetPrimaryVariableSwitches()
    {
        numPriVarsSwitched_ = 0;
        std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);
    }

    /*!
     * \brief Register all run-time parameters for the blackoil newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        BlackoilNewtonParams<Scalar>::registerParameters();
    }

    /*!
     * \brief Returns the number of degrees of freedom for which the
     *        interpretation has changed for the most recent iteration.
     */
    unsigned numPriVarsSwitched() const
    { return numPriVarsSwitched_; }

    const BlackoilNewtonParams<Scalar>& params() const { return bparams_; }
    const std::vector<std::uint8_t>& switchHistory() const { return wasSwitched_; }
    void setSwitchHistory(const std::vector<std::uint8_t>& history)
    {
        if (history.size() != wasSwitched_.size()) {
            OPM_THROW(std::invalid_argument, "Newton switch-history size mismatch");
        }
        wasSwitched_ = history;
    }
    void setNumPriVarsSwitched(unsigned count) { numPriVarsSwitched_ = count; }

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::beginIteration_
     */
    void beginIteration_()
    {
        numPriVarsSwitched_ = 0;
        ParentType::beginIteration_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     *
     * \param uCurrentIter Current solution iterator
     * \param uLastIter Last solution iterator
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
#if HAVE_MPI
        // in the MPI enabled case we need to add up the number of DOF
        // for which the interpretation changed over all processes.
        const int localSwitched = numPriVarsSwitched_;
        MPI_Allreduce(&localSwitched,
                      &numPriVarsSwitched_,
                      /*num=*/1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);
#endif // HAVE_MPI

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numPriVarsSwitched_;

        ParentType::endIteration_(uCurrentIter, uLastIter);
    }

public:
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& comm = this->simulator_.gridView().comm();

        int succeeded;
        try {
            ParentType::update_(nextSolution,
                                currentSolution,
                                solutionUpdate,
                                currentResidual);
            succeeded = 1;
        }
        catch (...) {
            succeeded = 0;
        }
        succeeded = comm.min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in adapting the primary variables");
        }

        numPriVarsSwitched_ = comm.sum(numPriVarsSwitched_);
    }

    template <class DofIndices>
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual,
                 const DofIndices& dofIndices)
    {
        const auto zero = 0.0 * solutionUpdate[0];
        for (auto dofIdx : dofIndices) {
            if (solutionUpdate[dofIdx] == zero) {
                continue;
            }
            updatePrimaryVariables_(dofIdx,
                                    nextSolution[dofIdx],
                                    currentSolution[dofIdx],
                                    solutionUpdate[dofIdx],
                                    currentResidual[dofIdx]);
        }
    }

protected:
    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        Valgrind::CheckDefined(currentResidual);
        wasSwitched_[globalDofIdx] = BlackOilNewtonUpdate<TypeTag>::update(
            this->problem(), FluidSystem{}, globalDofIdx, nextValue, currentValue,
            update, bparams_, wasSwitched_[globalDofIdx] != 0);
        numPriVarsSwitched_ += wasSwitched_[globalDofIdx] != 0;
    }

private:
    int numPriVarsSwitched_{};

    BlackoilNewtonParams<Scalar> bparams_{};

    // keep track of cells where the primary variable meaning has changed
    // to detect and hinder oscillations
    std::vector<std::uint8_t> wasSwitched_{};
};

} // namespace Opm

#endif // OPM_BLACK_OIL_NEWTHON_METHOD_HPP
