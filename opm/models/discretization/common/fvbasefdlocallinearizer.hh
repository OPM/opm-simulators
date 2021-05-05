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
 * \copydoc Opm::FvBaseFdLocalLinearizer
 */
#ifndef EWOMS_FV_BASE_FD_LOCAL_LINEARIZER_HH
#define EWOMS_FV_BASE_FD_LOCAL_LINEARIZER_HH

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Opm {
// forward declaration
template<class TypeTag>
class FvBaseFdLocalLinearizer;

} // namespace Opm

namespace Opm::Properties {

// declare the property tags required for the finite differences local linearizer

namespace TTag {
struct FiniteDifferenceLocalLinearizer {};
} // namespace TTag

template<class TypeTag, class MyTypeTag>
struct NumericDifferenceMethod { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct BaseEpsilon { using type = UndefinedProperty; };

// set the properties to be spliced in
template<class TypeTag>
struct LocalLinearizer<TypeTag, TTag::FiniteDifferenceLocalLinearizer>
{ using type = FvBaseFdLocalLinearizer<TypeTag>; };

template<class TypeTag>
struct Evaluation<TypeTag, TTag::FiniteDifferenceLocalLinearizer>
{ using type = GetPropType<TypeTag, Properties::Scalar>; };

/*!
 * \brief Specify which kind of method should be used to numerically
 * calculate the partial derivatives of the residual.
 *
 * -1 means backward differences, 0 means central differences, 1 means
 * forward differences. By default we use central differences.
 */
template<class TypeTag>
struct NumericDifferenceMethod<TypeTag, TTag::FiniteDifferenceLocalLinearizer> { static constexpr int value = +1; };

//! The base epsilon value for finite difference calculations
template<class TypeTag>
struct BaseEpsilon<TypeTag, TTag::FiniteDifferenceLocalLinearizer>
{
    using type = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr type value = std::max<type>(0.9123e-10, std::numeric_limits<type>::epsilon()*1.23e3);
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Calculates the Jacobian of the local residual for finite volume spatial
 *        discretizations using a finite difference method
 *
 * The local Jacobian for a given context is defined as the derivatives of the residuals
 * of all degrees of freedom featured by the stencil with regard to the primary variables
 * of the stencil's "primary" degrees of freedom.
 *
 * This class implements numeric differentiation using finite difference methods, i.e.
 * forward or backward differences (2nd order), or central differences (3rd order). The
 * method used is determined by the "NumericDifferenceMethod" property:
 *
 * - If the value of this property is smaller than 0, backward differences are used,
 *   i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - If the value of this property is 0, central differences are used, i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx
 *          \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward differences are used, i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx
 *          \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$ is the value of a
 * sub-control volume's primary variable at the evaluation point and \f$\epsilon\f$ is a
 * small scalar value larger than 0.
 */
template<class TypeTag>
class FvBaseFdLocalLinearizer
{
private:
    using Implementation = GetPropType<TypeTag, Properties::LocalLinearizer>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };

    // extract local matrices from jacobian matrix for consistency
    using ScalarMatrixBlock = typename GetPropType<TypeTag, Properties::SparseMatrixAdapter>::MatrixBlock;
    using ScalarVectorBlock = Dune::FieldVector<Scalar, numEq>;

    using ScalarLocalBlockVector = Dune::BlockVector<ScalarVectorBlock>;
    using ScalarLocalBlockMatrix = Dune::Matrix<ScalarMatrixBlock>;

    using LocalEvalBlockVector = typename LocalResidual::LocalEvalBlockVector;

#if __GNUC__ == 4 && __GNUC_MINOR__ <= 6
public:
    // make older GCCs happy by providing a public copy constructor (this is necessary
    // for their implementation of std::vector, although the method is never called...)
    FvBaseFdLocalLinearizer(const FvBaseFdLocalLinearizer&)
        : internalElemContext_(0)
    {}

#else
    // copying local residual objects around is a very bad idea, so we explicitly prevent
    // it...
    FvBaseFdLocalLinearizer(const FvBaseFdLocalLinearizer&) = delete;
#endif
public:
    FvBaseFdLocalLinearizer()
        : internalElemContext_(0)
    { }

    ~FvBaseFdLocalLinearizer()
    { delete internalElemContext_; }

    /*!
     * \brief Register all run-time parameters for the local jacobian.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, NumericDifferenceMethod,
                             "The method used for numeric differentiation (-1: backward "
                             "differences, 0: central differences, 1: forward differences)");
    }

    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param simulator The simulator object of the simulation.
     */
    void init(Simulator& simulator)
    {
        simulatorPtr_ = &simulator;
        delete internalElemContext_;
        internalElemContext_ = new ElementContext(simulator);
    }

    /*!
     * \brief Compute an element's local Jacobian matrix and evaluate its residual.
     *
     * The local Jacobian for a given context is defined as the derivatives of the
     * residuals of all degrees of freedom featured by the stencil with regard to the
     * primary variables of the stencil's "primary" degrees of freedom. Adding the local
     * Jacobians for all elements in the grid will give the global Jacobian 'grad f(x)'.
     *
     * \param element The grid element for which the local residual and its local
     *                Jacobian should be calculated.
     */
    void linearize(const Element& element)
    {
        linearize(*internalElemContext_, element);
    }

    /*!
     * \brief Compute an element's local Jacobian matrix and evaluate its residual.
     *
     * The local Jacobian for a given context is defined as the derivatives of the
     * residuals of all degrees of freedom featured by the stencil with regard to the
     * primary variables of the stencil's "primary" degrees of freedom. Adding the local
     * Jacobians for all elements in the grid will give the global Jacobian 'grad f(x)'.
     *
     * After calling this method the ElementContext is in an undefined state, so do not
     * use it anymore!
     *
     * \param elemCtx The element execution context for which the local residual and its
     *                local Jacobian should be calculated.
     */
    void linearize(ElementContext& elemCtx, const Element& elem)
    {
        elemCtx.updateAll(elem);

        // update the weights of the primary variables for the context
        model_().updatePVWeights(elemCtx);

        resize_(elemCtx);
        reset_(elemCtx);

        // calculate the local residual
        localResidual_.eval(residual_, elemCtx);

        // calculate the local jacobian matrix
        size_t numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        for (unsigned dofIdx = 0; dofIdx < numPrimaryDof; dofIdx++) {
            for (unsigned pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().evalPartialDerivative_(elemCtx, dofIdx, pvIdx);

                // incorporate the partial derivatives into the local Jacobian matrix
                updateLocalJacobian_(elemCtx, dofIdx, pvIdx);
            }
        }
    }

    /*!
     * \brief Returns the unweighted epsilon value used to calculate
     *        the local derivatives
     */
    static Scalar baseEpsilon()
    { return getPropValue<TypeTag, Properties::BaseEpsilon>(); }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemCtx The element execution context for which the
     *                local residual and its gradient should be
     *                calculated.
     * \param dofIdx     The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon(const ElementContext& elemCtx,
                          unsigned dofIdx,
                          unsigned pvIdx) const
    {
        unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        Scalar pvWeight = elemCtx.model().primaryVarWeight(globalIdx, pvIdx);
        assert(pvWeight > 0 && std::isfinite(pvWeight));
        Valgrind::CheckDefined(pvWeight);

        return baseEpsilon()/pvWeight;
    }

    /*!
     * \brief Return reference to the local residual.
     */
    LocalResidual& localResidual()
    { return localResidual_; }

    /*!
     * \brief Return reference to the local residual.
     */
    const LocalResidual& localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns the local Jacobian matrix of the residual of a sub-control volume.
     *
     * \param domainScvIdx The local index of the sub control volume
     *                     which contains the independents
     * \param rangeScvIdx The local index of the sub control volume
     *                    which contains the local residual
     */
    const ScalarMatrixBlock& jacobian(unsigned domainScvIdx, unsigned rangeScvIdx) const
    { return jacobian_[domainScvIdx][rangeScvIdx]; }

    /*!
     * \brief Returns the local residual of a sub-control volume.
     *
     * \param dofIdx The local index of the sub control volume
     */
    const ScalarVectorBlock& residual(unsigned dofIdx) const
    { return residual_[dofIdx]; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Simulator& simulator_() const
    { return *simulatorPtr_; }
    const Problem& problem_() const
    { return simulatorPtr_->problem(); }
    const Model& model_() const
    { return simulatorPtr_->model(); }

    /*!
     * \brief Returns the numeric difference method which is applied.
     */
    static int numericDifferenceMethod_()
    { return EWOMS_GET_PARAM(TypeTag, int, NumericDifferenceMethod); }

    /*!
     * \brief Resize all internal attributes to the size of the
     *        element.
     */
    void resize_(const ElementContext& elemCtx)
    {
        size_t numDof = elemCtx.numDof(/*timeIdx=*/0);
        size_t numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);

        residual_.resize(numDof);
        jacobian_.setSize(numDof, numPrimaryDof);

        derivResidual_.resize(numDof);
    }

    /*!
     * \brief Reset the all relevant internal attributes to 0
     */
    void reset_(const ElementContext& elemCtx)
    {
        size_t numDof = elemCtx.numDof(/*timeIdx=*/0);
        size_t numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        for (unsigned primaryDofIdx = 0; primaryDofIdx < numPrimaryDof; ++ primaryDofIdx)
            for (unsigned dof2Idx = 0; dof2Idx < numDof; ++ dof2Idx)
                jacobian_[dof2Idx][primaryDofIdx] = 0.0;

        for (unsigned primaryDofIdx = 0; primaryDofIdx < numDof; ++ primaryDofIdx)
            residual_[primaryDofIdx] = 0.0;
    }

    /*!
     * \brief Compute the partial derivatives of a context's residual functions
     *
     * This method can be overwritten by the implementation if a better scheme than
     * numerical differentiation is available.
     *
     * The default implementation of this method uses numeric differentiation,
     * i.e. forward or backward differences (2nd order), or central differences (3rd
     * order). The method used is determined by the "NumericDifferenceMethod" property:
     *
     * - If the value of this property is smaller than 0, backward differences are used,
     *   i.e.:
     *   \f[
     *     \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
     *   \f]
     *
     * - If the value of this property is 0, central differences are used, i.e.:
     *   \f[
     *     \frac{\partial f(x)}{\partial x} \approx
     *          \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
     *   \f]
     *
     * - if the value of this property is larger than 0, forward
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
     *   \f]
     *
     * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$ is the value
     * of a sub-control volume's primary variable at the evaluation point and
     * \f$\epsilon\f$ is a small value larger than 0.
     *
     * \param elemCtx The element context for which the local partial
     *                derivative ought to be calculated
     * \param dofIdx The sub-control volume index of the current
     *               finite element for which the partial derivative
     *               ought to be calculated
     * \param pvIdx The index of the primary variable at the dofIdx'
     *              sub-control volume of the current finite element
     *              for which the partial derivative ought to be
     *              calculated
     */
    void evalPartialDerivative_(ElementContext& elemCtx,
                                unsigned dofIdx,
                                unsigned pvIdx)
    {
        // save all quantities which depend on the specified primary
        // variable at the given sub control volume
        elemCtx.stashIntensiveQuantities(dofIdx);

        PrimaryVariables priVars(elemCtx.primaryVars(dofIdx, /*timeIdx=*/0));
        Scalar eps = asImp_().numericEpsilon(elemCtx, dofIdx, pvIdx);
        Scalar delta = 0.0;

        if (numericDifferenceMethod_() >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the deflected residual
            elemCtx.updateIntensiveQuantities(priVars, dofIdx, /*timeIdx=*/0);
            elemCtx.updateAllExtensiveQuantities();
            localResidual_.eval(derivResidual_, elemCtx);
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ = residual_;
        }

        if (numericDifferenceMethod_() <= 0) {
            // we are not using forward differences, i.e. we don't
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // calculate the deflected residual again, this time we use the local
            // residual's internal storage.
            elemCtx.updateIntensiveQuantities(priVars, dofIdx, /*timeIdx=*/0);
            elemCtx.updateAllExtensiveQuantities();
            localResidual_.eval(elemCtx);

            derivResidual_ -= localResidual_.residual();
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ -= residual_;
        }

        assert(delta > 0);

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        derivResidual_ /= delta;

        // restore the original state of the element's volume
        // variables
        elemCtx.restoreIntensiveQuantities(dofIdx);

#ifndef NDEBUG
        for (unsigned i = 0; i < derivResidual_.size(); ++i)
            Valgrind::CheckDefined(derivResidual_[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the partial derivatives of
     *        all equations for primary variable 'pvIdx' at the degree of freedom
     *        associated with 'focusDofIdx'.
     */
    void updateLocalJacobian_(const ElementContext& elemCtx,
                              unsigned focusDofIdx,
                              unsigned pvIdx)
    {
        size_t numDof = elemCtx.numDof(/*timeIdx=*/0);
        for (unsigned dofIdx = 0; dofIdx < numDof; dofIdx++) {
            for (unsigned eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[dofIdx][focusDofIdx][eqIdx][pvIdx] is the partial derivative of the
                // residual function 'eqIdx' for the degree of freedom 'dofIdx' with
                // regard to the primary variable 'pvIdx' of the degree of freedom
                // 'focusDofIdx'
                jacobian_[dofIdx][focusDofIdx][eqIdx][pvIdx] = derivResidual_[dofIdx][eqIdx];
                Valgrind::CheckDefined(jacobian_[dofIdx][focusDofIdx][eqIdx][pvIdx]);
            }
        }
    }

    Simulator *simulatorPtr_;
    Model *modelPtr_;

    ElementContext *internalElemContext_;

    LocalEvalBlockVector residual_;
    LocalEvalBlockVector derivResidual_;
    ScalarLocalBlockMatrix jacobian_;

    LocalResidual localResidual_;
};

} // namespace Opm

#endif
