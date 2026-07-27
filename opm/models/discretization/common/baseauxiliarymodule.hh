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
#include <opm/simulators/linalg/linalgproperties.hh>

#include <set>
#include <vector>

namespace Opm::Properties::Tag {

struct AuxModule {};

} // namespace Opm::Properties::TTag

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
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

protected:
    using NeighborSet = std::set<unsigned>;

public:
    /*!
     * \brief A flux connection between an auxiliary degree of freedom and another
     *        degree of freedom of the model.
     *
     * This is for auxiliary modules whose degrees of freedom carry the *model's own*
     * conservation equations -- an auxiliary "cell" with an authored volume and an
     * authored connection list, as opposed to a genuinely different unknown such as a
     * well's bottom hole pressure. Such a module only has to declare which pairs of
     * degrees of freedom exchange fluxes; the discretization then assembles them with
     * the same local residual it uses for the grid, reading the transmissibility (and
     * the thermal/diffusive counterparts) from the problem exactly as it does for a
     * geometric face. Both endpoints are stored as plain degree-of-freedom indices, so
     * a connection may join two auxiliary degrees of freedom or an auxiliary one and a
     * grid cell.
     */
    struct AuxiliaryConnection
    {
        unsigned dof1{};
        unsigned dof2{};
    };

    virtual ~BaseAuxiliaryModule() = default;

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
    int dofOffset() const
    { return dofOffset_; }

    /*!
     * \brief Whether this module's degrees of freedom carry the model's own
     *        conservation equations.
     *
     * False by default, which describes a module whose unknowns are of a different
     * kind -- a well's bottom hole pressure, a mortar multiplier -- and which
     * assembles and scales its own equations in linearize().  Such degrees of freedom
     * are deliberately kept out of the model's error norm and primary-variable
     * switching.
     *
     * A module which returns true is declaring the opposite: its degrees of freedom are
     * cells as far as the model is concerned, with the model's own unknowns, and they
     * take part in the Newton update, the convergence measures and the variable
     * switching exactly like a grid cell.
     */
    virtual bool carriesModelEquations() const
    { return false; }

    /*!
     * \brief The volume associated with one of this module's degrees of freedom.
     *
     * Zero unless the module's degrees of freedom are cells.  A grid degree of freedom
     * takes this from the geometry of its entity; an auxiliary cell has no entity, so it
     * has to state the volume itself.
     */
    virtual Scalar dofVolume(unsigned /*localDofIdx*/) const
    { return 0.0; }

    /*!
     * \brief Given a degree of freedom relative to the current auxiliary equation,
     *        return the corresponding index in the global system of equations.
     */
    int localToGlobalDof(unsigned localDofIdx) const
    {
        assert(localDofIdx < numDofs());
        return dofOffset_ + localDofIdx;
    }

    /*!
     * \brief Specify the additional neighboring correlations caused by the auxiliary
     *        module.
     */
    virtual void addNeighbors(std::vector<NeighborSet>& neighbors) const = 0;

    /*!
     * \brief Append this module's flux connections, if any.
     *
     * Modules whose degrees of freedom do not carry the model's conservation equations
     * -- the well models, for instance, which assemble their own equations in
     * linearize() -- leave this empty, which is the default.
     *
     * The discretization inserts each reported connection into the sparsity pattern in
     * both directions and assembles it from both endpoints, so a connection must be
     * reported exactly once, not once per endpoint.
     */
    virtual void addConnections(std::vector<AuxiliaryConnection>&) const
    {}

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
    virtual void postSolve(GlobalEqVector&)
    {}

private:
    int dofOffset_{};
};

} // namespace Opm

#endif
