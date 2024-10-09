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
 * \copydoc Opm::FlowExpNewtonMethod
 */
#ifndef OPM_FLOW_EXP_NEWTON_METHOD_HPP
#define OPM_FLOW_EXP_NEWTON_METHOD_HPP

#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/blackoil/blackoilnewtonmethod.hpp>
#include <opm/models/utils/signum.hh>

namespace Opm::Parameters {

// the tolerated amount of "incorrect" amount of oil per time step for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
template<class Scalar>
struct EclNewtonSumTolerance { static constexpr Scalar value = 1e-5; };

// make all Newton iterations strict, i.e., the volumetric Newton tolerance must be
// always be upheld in the majority of the spatial domain. In this context, "majority"
// means 1 - EclNewtonRelaxedVolumeFraction.
struct EclNewtonStrictIterations { static constexpr int value = 100; };

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
template<class Scalar>
struct EclNewtonRelaxedVolumeFraction { static constexpr Scalar value = 0.05; };

template<class Scalar>
struct EclNewtonSumToleranceExponent { static constexpr Scalar value = 1.0 / 3.0; };

// the maximum volumetric error of a cell in the relaxed region
template<class Scalar>
struct EclNewtonRelaxedTolerance { static constexpr Scalar value = NewtonTolerance<Scalar>::value * 1e6; };

} // namespace Opm::Parameters

namespace Opm {

/*!
 * \brief A newton solver.
 */
template <class TypeTag>
class FlowExpNewtonMethod : public BlackOilNewtonMethod<TypeTag>
{
    using ParentType = BlackOilNewtonMethod<TypeTag>;
    using DiscNewtonMethod = GetPropType<TypeTag, Properties::DiscNewtonMethod>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

    static constexpr int contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr int contiEnergyEqIdx = Indices::contiEnergyEqIdx;

    friend NewtonMethod<TypeTag>;
    friend DiscNewtonMethod;
    friend ParentType;

public:
    explicit FlowExpNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        errorPvFraction_ = 1.0;
        relaxedMaxPvFraction_ = Parameters::Get<Parameters::EclNewtonRelaxedVolumeFraction<Scalar>>();

        sumTolerance_ = 0.0; // this gets determined in the error calculation proceedure
        relaxedTolerance_ = Parameters::Get<Parameters::EclNewtonRelaxedTolerance<Scalar>>();

        numStrictIterations_ = Parameters::Get<Parameters::EclNewtonStrictIterations>();
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        Parameters::Register<Parameters::EclNewtonSumTolerance<Scalar>>
                ("The maximum error tolerated by the Newton "
                 "method for considering a solution to be converged");
        Parameters::Register<Parameters::EclNewtonStrictIterations>
                 ("The number of Newton iterations where the "
                  "volumetric error is considered.");
        Parameters::Register<Parameters::EclNewtonRelaxedVolumeFraction<Scalar>>
                  ("The fraction of the pore volume of the reservoir "
                   "where the volumetric error may be violated during strict Newton iterations.");
        Parameters::Register<Parameters::EclNewtonSumToleranceExponent<Scalar>>
                  ("The the exponent used to scale the sum tolerance by "
                   "the total pore volume of the reservoir.");
        Parameters::Register<Parameters::EclNewtonRelaxedTolerance<Scalar>>
                   ("The maximum error which the volumetric residual "
                     "may exhibit if it is in a 'relaxed' region during a strict iteration.");
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    {
        if (errorPvFraction_ < relaxedMaxPvFraction_) {
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;
        } else if (this->numIterations() > numStrictIterations_) {
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;
        }

        return this->error_ <= this->tolerance() && errorSum_ <= sumTolerance_;
    }

    void preSolve_(const SolutionVector&,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = this->model().linearizer().constraintsMap();
        this->lastError_ = this->error_;
        Scalar newtonMaxError = this->params_.maxError_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        this->error_ = 0.0;
        Dune::FieldVector<Scalar, numEq> componentSumError;
        std::fill(componentSumError.begin(), componentSumError.end(), 0.0);
        Scalar sumPv = 0.0;
        errorPvFraction_ = 0.0;
        const Scalar dt = this->simulator_.timeStepSize();
        for (unsigned dofIdx = 0; dofIdx < currentResidual.size(); ++dofIdx) {
            // do not consider auxiliary DOFs for the error
            if (dofIdx >= this->model().numGridDof()
                || this->model().dofTotalVolume(dofIdx) <= 0.0) {
                continue;
            }

            if (!this->model().isLocalDof(dofIdx)) {
                continue;
            }

            // also do not consider DOFs which are constraint
            if (this->enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0) {
                    continue;
                }
            }

            const auto& r = currentResidual[dofIdx];
            Scalar pvValue =   this->simulator_.problem().referencePorosity(dofIdx, /*timeIdx=*/0) *
                               this->model().dofTotalVolume(dofIdx);
            sumPv += pvValue;
            bool cnvViolated = false;

            Scalar dofVolume = this->model().dofTotalVolume(dofIdx);

            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                Scalar tmpError = r[eqIdx] * dt * this->model().eqWeight(dofIdx, eqIdx) / pvValue;
                Scalar tmpError2 = r[eqIdx] * this->model().eqWeight(dofIdx, eqIdx);

                // in the case of a volumetric formulation, the residual in the above is
                // per cubic meter
                if (getPropValue<TypeTag, Properties::UseVolumetricResidual>()) {
                    tmpError *= dofVolume;
                    tmpError2 *= dofVolume;
                }

                this->error_ = max(std::abs(tmpError), this->error_);

                if (std::abs(tmpError) > this->params_.tolerance_) {
                    cnvViolated = true;
                }

                componentSumError[eqIdx] += std::abs(tmpError2);
            }
            if (cnvViolated) {
                errorPvFraction_ += pvValue;
            }
        }

        // take the other processes into account
        this->error_ = this->comm_.max(this->error_);
        componentSumError = this->comm_.sum(componentSumError);
        sumPv = this->comm_.sum(sumPv);
        errorPvFraction_ = this->comm_.sum(errorPvFraction_);

        componentSumError /= sumPv;
        componentSumError *= dt;

        errorPvFraction_ /= sumPv;

        errorSum_ = 0;
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            errorSum_ = std::max(std::abs(componentSumError[eqIdx]), errorSum_);
        }

        // scale the tolerance for the total error with the pore volume. by default, the
        // exponent is 1/3, i.e., cubic root.
        Scalar x = Parameters::Get<Parameters::EclNewtonSumTolerance<Scalar>>();
        Scalar y = Parameters::Get<Parameters::EclNewtonSumToleranceExponent<Scalar>>();
        sumTolerance_ = x*std::pow(sumPv, y);

        this->endIterMsg() << " (max: " << this->params_.tolerance_
                           << ", violated for " << errorPvFraction_ * 100
                           << "% of the pore volume), aggegate error: "
                           << errorSum_ << " (max: " << sumTolerance_ << ")";

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (this->error_ > newtonMaxError) {
            throw NumericalProblem("Newton: Error "+std::to_string(double(this->error_))
                                   + " is larger than maximum allowed error of "
                                   + std::to_string(double(newtonMaxError)));
        }

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (errorSum_ > newtonMaxError) {
            throw NumericalProblem("Newton: Sum of the error "+std::to_string(double(errorSum_))
                                   + " is larger than maximum allowed error of "
                                   + std::to_string(double(newtonMaxError)));
        }
    }

    void endIteration_(SolutionVector& nextSolution,
                       const SolutionVector& currentSolution)
    {
        ParentType::endIteration_(nextSolution, currentSolution);
        OpmLog::debug( "Newton iteration " + std::to_string(this->numIterations_) + ""
                  + " error: " + std::to_string(double(this->error_))
                  + this->endIterMsg().str());
        this->endIterMsg().str("");
    }

private:
    Scalar errorPvFraction_;
    Scalar errorSum_;

    Scalar relaxedTolerance_;
    Scalar relaxedMaxPvFraction_;

    Scalar sumTolerance_;

    int numStrictIterations_;
};

} // namespace Opm

#endif
