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
 * \copydoc Opm::FvBaseNewtonMethod
 */
#ifndef EWOMS_FV_BASE_NEWTON_METHOD_HH
#define EWOMS_FV_BASE_NEWTON_METHOD_HH

#include "fvbasenewtonconvergencewriter.hh"

#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

template <class TypeTag>
class FvBaseNewtonMethod;

template <class TypeTag>
class FvBaseNewtonConvergenceWriter;
} // namespace Opm

namespace Opm::Properties {

//! create a type tag for the Newton method of the finite-volume discretization
// Create new type tags
namespace TTag {
struct FvBaseNewtonMethod { using InheritsFrom = std::tuple<NewtonMethod>; };
} // end namespace TTag

//! The discretization specific part of he implementing the Newton algorithm
template<class TypeTag, class MyTypeTag>
struct DiscNewtonMethod { using type = UndefinedProperty; };

// set default values
template<class TypeTag>
struct DiscNewtonMethod<TypeTag, TTag::FvBaseNewtonMethod>
{ using type = FvBaseNewtonMethod<TypeTag>; };

template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::FvBaseNewtonMethod>
{ using type = GetPropType<TypeTag, Properties::DiscNewtonMethod>; };

template<class TypeTag>
struct NewtonConvergenceWriter<TypeTag, TTag::FvBaseNewtonMethod>
{ using type = FvBaseNewtonConvergenceWriter<TypeTag>; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief A Newton method for models using a finite volume discretization.
 *
 * This class is sufficient for most models which use an Element or a
 * Vertex Centered Finite Volume discretization.
 */
template <class TypeTag>
class FvBaseNewtonMethod : public NewtonMethod<TypeTag>
{
    using ParentType = NewtonMethod<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::NewtonMethod>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;


public:
    FvBaseNewtonMethod(Simulator& simulator)
        : ParentType(simulator)
    { }

protected:
    friend class NewtonMethod<TypeTag>;

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the converged() and
     * proceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param nextSolution The solution vector at the end of the current iteration
     * \param currentSolution The solution vector at the beginning of the current iteration
     * \param solutionUpdate The delta as calculated by solving the linear system of
     *                       equations. This parameter also stores the updated solution.
     * \param currentResidual The residual (i.e., right-hand-side) of the current solution.
     */
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        ParentType::update_(nextSolution, currentSolution, solutionUpdate, currentResidual);

        // make sure that the intensive quantities get recalculated at the next
        // linearization
        if (model_().storeIntensiveQuantities()) {
            for (unsigned dofIdx = 0; dofIdx < model_().numGridDof(); ++dofIdx)
                model_().setIntensiveQuantitiesCacheEntryValidity(dofIdx,
                                                                  /*timeIdx=*/0,
                                                                  /*valid=*/false);
        }
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void beginIteration_()
    {
        model_().syncOverlap();

        ParentType::beginIteration_();
    }

    /*!
     * \brief Returns a reference to the model.
     */
    Model& model_()
    { return ParentType::model(); }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model& model_() const
    { return ParentType::model(); }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Opm

#endif
