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
 * \copydoc Opm::EclNewtonMethod
 */
#ifndef EWOMS_ECL_NEWTON_METHOD_SEQ_HH
#define EWOMS_ECL_NEWTON_METHOD_SEQ_HH

#include <ebos/eclblackoilnewtonmethod.hh>
#include <opm/models/utils/signum.hh>
#include <opm/common/OpmLog/OpmLog.hpp>


#include <opm/material/common/Unused.hpp>

namespace Opm {

/*!
 * \brief A newton solver which is ebos specific.
 */
template <class TypeTag>
class EclNewtonMethodSeq : public EclBlackOilNewtonMethod<TypeTag>
{
    typedef EclBlackOilNewtonMethod<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) DiscNewtonMethod;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    static const unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

    static constexpr int contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr int contiEnergyEqIdx = Indices::contiEnergyEqIdx;

    friend NewtonMethod<TypeTag>;
    friend DiscNewtonMethod;
    friend ParentType;

public:
    EclNewtonMethodSeq(Simulator& simulator) : ParentType(simulator)
    {
        errorPvFraction_ = 1.0;
        relaxedMaxPvFraction_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonRelaxedVolumeFraction);

        sumTolerance_ = 0.0; // this gets determined in the error calculation proceedure
        relaxedTolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonRelaxedTolerance);

        numStrictIterations_ = EWOMS_GET_PARAM(TypeTag, int, EclNewtonStrictIterations);
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonSumTolerance,
                             "The maximum error tolerated by the Newton"
                             "method for considering a solution to be "
                             "converged");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclNewtonStrictIterations,
                             "The number of Newton iterations where the"
                             " volumetric error is considered.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonRelaxedVolumeFraction,
                             "The fraction of the pore volume of the reservoir "
                             "where the volumetric error may be voilated during "
                             "strict Newton iterations.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonSumToleranceExponent,
                             "The the exponent used to scale the sum tolerance by "
                             "the total pore volume of the reservoir.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonRelaxedTolerance,
                             "The maximum error which the volumetric residual "
                             "may exhibit if it is in a 'relaxed' "
                             "region during a strict iteration.");
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    {
        if (errorPvFraction_ < relaxedMaxPvFraction_)
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;
        else if (this->numIterations() > numStrictIterations_)
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;

        return this->error_ <= this->tolerance() && errorSum_ <= sumTolerance_;
    }

    void preSolve_(const SolutionVector& currentSolution  OPM_UNUSED,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = this->model().linearizer().constraintsMap();
        this->lastError_ = this->error_;
        Scalar newtonMaxError = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError);

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
                || this->model().dofTotalVolume(dofIdx) <= 0.0)
                continue;

            if (!this->model().isLocalDof(dofIdx))
                continue;

            // also do not consider DOFs which are constraint
            if (this->enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0)
                    continue;
            }

            const auto& r = currentResidual[dofIdx];
            Scalar pvValue =
                this->simulator_.problem().referencePorosity(dofIdx, /*timeIdx=*/0)
                * this->model().dofTotalVolume(dofIdx);
            sumPv += pvValue;
            bool cnvViolated = false;

            Scalar dofVolume = this->model().dofTotalVolume(dofIdx);

            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                Scalar tmpError = r[eqIdx] * dt * this->model().eqWeight(dofIdx, eqIdx) / pvValue;
                Scalar tmpError2 = r[eqIdx] * this->model().eqWeight(dofIdx, eqIdx);

                // in the case of a volumetric formulation, the residual in the above is
                // per cubic meter
                if (GET_PROP_VALUE(TypeTag, UseVolumetricResidual)) {
                    tmpError *= dofVolume;
                    tmpError2 *= dofVolume;
                }

                this->error_ = Opm::max(std::abs(tmpError), this->error_);

                if (std::abs(tmpError) > this->tolerance_)
                    cnvViolated = true;

                componentSumError[eqIdx] += std::abs(tmpError2);
            }
            if (cnvViolated)
                errorPvFraction_ += pvValue;
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
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
            errorSum_ = std::max(std::abs(componentSumError[eqIdx]), errorSum_);

        // scale the tolerance for the total error with the pore volume. by default, the
        // exponent is 1/3, i.e., cubic root.
        Scalar x = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonSumTolerance);
        Scalar y = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonSumToleranceExponent);
        sumTolerance_ = x*std::pow(sumPv, y);

        this->endIterMsg() << " (max: " << this->tolerance_ << ", violated for " << errorPvFraction_*100 << "% of the pore volume), aggegate error: " << errorSum_ << " (max: " << sumTolerance_ << ")";

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (this->error_ > newtonMaxError)
            throw Opm::NumericalIssue("Newton: Error "+std::to_string(double(this->error_))
                                        +" is larger than maximum allowed error of "
                                        +std::to_string(double(newtonMaxError)));

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (errorSum_ > newtonMaxError)
            throw Opm::NumericalIssue("Newton: Sum of the error "+std::to_string(double(errorSum_))
                                        +" is larger than maximum allowed error of "
                                        +std::to_string(double(newtonMaxError)));
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
