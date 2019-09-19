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
 * \copydoc Opm::BaseAuxiliaryModule
 */
#ifndef EWOMS_BASE_AUXILIARY_MODULE_HH
#define EWOMS_BASE_AUXILIARY_MODULE_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <set>
#include <vector>

BEGIN_PROPERTIES

NEW_TYPE_TAG(AuxModule);

// declare the properties required by the for the ecl grid manager
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(DofMapper);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(SparseMatrixAdapter);

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup ModelModules
 *
 * \brief Base class for specifying auxiliary equations.
 *
 * For example, these equations can be wells, non-neighboring connections, interfaces
 * between model domains, etc.
 */
template <class TypeTag>
class BaseAuxiliaryModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

protected:
    typedef std::set<unsigned> NeighborSet;

public:
    virtual ~BaseAuxiliaryModule()
    {}

    /*!
     * \brief Returns the number of additional degrees of freedom required for the
     *        auxiliary module.
     */
    virtual unsigned numDofs() const = 0;

    /*!
     * \brief Set the offset in the global system of equations for the first degree of
     *        freedom of this auxiliary module.
     */
    void setDofOffset(int value)
    { dofOffset_ = value; }

    /*!
     * \brief Return the offset in the global system of equations for the first degree of
     *        freedom of this auxiliary module.
     */
    int dofOffset()
    { return dofOffset_; }

    /*!
     * \brief Given a degree of freedom relative to the current auxiliary equation,
     *        return the corresponding index in the global system of equations.
     */
    int localToGlobalDof(unsigned localDofIdx) const
    {
        assert(0 <= localDofIdx);
        assert(localDofIdx < numDofs());
        return dofOffset_ + localDofIdx;
    }

    /*!
     * \brief Specify the additional neighboring correlations caused by the auxiliary
     *        module.
     */
    virtual void addNeighbors(std::vector<NeighborSet>& neighbors) const = 0;

    /*!
     * \brief Set the initial condition of the auxiliary module in the solution vector.
     */
    virtual void applyInitial() = 0;

    /*!
     * \brief Linearize the auxiliary equation.
     */
    virtual void linearize(SparseMatrixAdapter& matrix, GlobalEqVector& residual) = 0;

    /*!
     * \brief This method is called after the linear solver has been called but before
     *        the solution is updated for the next iteration.
     *
     * It is intended to implement stuff like Schur complements.
     */
    virtual void postSolve(GlobalEqVector& residual OPM_UNUSED)
    {};

private:
    int dofOffset_;
};

} // namespace Opm

#endif
