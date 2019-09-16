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
 * \brief Provides data handles for parallel communication which
 *        operate on DOFs
 */
#ifndef EWOMS_GRID_COMM_HANDLES_HH
#define EWOMS_GRID_COMM_HANDLES_HH

#include <opm/material/common/Unused.hpp>

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/version.hh>

namespace Opm {

/*!
 * \brief Data handle for parallel communication which sums up all
 *        values are attached to DOFs
 */
template <class FieldType, class Container, class EntityMapper, int commCodim>
class GridCommHandleSum
    : public Dune::CommDataHandleIF<GridCommHandleSum<FieldType, Container,
                                                      EntityMapper, commCodim>,
                                    FieldType>
{
public:
    GridCommHandleSum(Container& container, const EntityMapper& mapper)
        : mapper_(mapper), container_(container)
    {}

    bool contains(int dim OPM_UNUSED, int codim) const
    {
        // return true if the codim is the same as the codim which we
        // are asked to communicate with.
        return codim == commCodim;
    }

    bool fixedsize(int dim OPM_UNUSED, int codim OPM_UNUSED) const
    {
        // for each DOF we communicate a single value which has a
        // fixed size
        return true;
    }

    template <class EntityType>
    size_t size(const EntityType& e OPM_UNUSED) const
    {
        // communicate a field type per entity
        return 1;
    }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.write(container_[dofIdx]);
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n OPM_UNUSED)
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));

        FieldType tmp;
        buff.read(tmp);
        container_[dofIdx] += tmp;
    }

private:
    const EntityMapper& mapper_;
    Container& container_;
};

/*!
 * \brief Data handle for parallel communication which can be used to
 *        set the values values of ghost and overlap DOFs from their
 *        respective master processes.
 */
template <class FieldType, class Container, class EntityMapper, unsigned commCodim>
class GridCommHandleGhostSync
    : public Dune::CommDataHandleIF<GridCommHandleGhostSync<FieldType, Container,
                                                            EntityMapper, commCodim>,
                                    FieldType>
{
public:
    GridCommHandleGhostSync(Container& container, const EntityMapper& mapper)
        : mapper_(mapper), container_(container)
    {
    }

    bool contains(unsigned dim OPM_UNUSED, unsigned codim) const
    {
        // return true if the codim is the same as the codim which we
        // are asked to communicate with.
        return codim == commCodim;
    }

    bool fixedsize(unsigned dim OPM_UNUSED, unsigned codim OPM_UNUSED) const
    {
        // for each DOF we communicate a single value which has a
        // fixed size
        return true;
    }

    template <class EntityType>
    size_t size(const EntityType& e OPM_UNUSED) const
    {
        // communicate a field type per entity
        return 1;
    }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.write(container_[dofIdx]);
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n OPM_UNUSED)
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.read(container_[dofIdx]);
    }

private:
    const EntityMapper& mapper_;
    Container& container_;
};

/*!
 * \brief Data handle for parallel communication which takes the
 *        maximum of all values that are attached to DOFs
 */
template <class FieldType, class Container, class EntityMapper, unsigned commCodim>
class GridCommHandleMax
    : public Dune::CommDataHandleIF<GridCommHandleMax<FieldType, Container,
                                                      EntityMapper, commCodim>,
                                    FieldType>
{
public:
    GridCommHandleMax(Container& container, const EntityMapper& mapper)
        : mapper_(mapper), container_(container)
    {}

    bool contains(unsigned dim OPM_UNUSED, unsigned codim) const
    {
        // return true if the codim is the same as the codim which we
        // are asked to communicate with.
        return codim == commCodim;
    }

    bool fixedsize(unsigned dim OPM_UNUSED, unsigned codim OPM_UNUSED) const
    {
        // for each DOF we communicate a single value which has a
        // fixed size
        return true;
    }

    template <class EntityType>
    size_t size(const EntityType& e OPM_UNUSED) const
    {
        // communicate a field type per entity
        return 1;
    }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.write(container_[dofIdx]);
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n OPM_UNUSED)
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        FieldType tmp;
        buff.read(tmp);
        container_[dofIdx] = std::max(container_[dofIdx], tmp);
    }

private:
    const EntityMapper& mapper_;
    Container& container_;
};

/*!
 * \brief Provides data handle for parallel communication which takes
 *        the minimum of all values that are attached to DOFs
 */
template <class FieldType, class Container, class EntityMapper, unsigned commCodim>
class GridCommHandleMin
    : public Dune::CommDataHandleIF<GridCommHandleMin<FieldType, Container,
                                                      EntityMapper, commCodim>,
                                    FieldType>
{
public:
    GridCommHandleMin(Container& container, const EntityMapper& mapper)
        : mapper_(mapper), container_(container)
    {}

    bool contains(unsigned dim OPM_UNUSED, unsigned codim) const
    {
        // return true if the codim is the same as the codim which we
        // are asked to communicate with.
        return codim == commCodim;
    }

    bool fixedsize(unsigned dim OPM_UNUSED, unsigned codim OPM_UNUSED) const
    {
        // for each DOF we communicate a single value which has a
        // fixed size
        return true;
    }

    template <class EntityType>
    size_t size(const EntityType& e OPM_UNUSED) const
    {
        // communicate a field type per entity
        return 1;
    }

    template <class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp& buff, const EntityType& e) const
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        buff.write(container_[dofIdx]);
    }

    template <class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp& buff, const EntityType& e, size_t n OPM_UNUSED)
    {
        unsigned dofIdx = static_cast<unsigned>(mapper_.index(e));
        FieldType tmp;
        buff.read(tmp);
        container_[dofIdx] = std::min(container_[dofIdx], tmp);
    }

private:
    const EntityMapper& mapper_;
    Container& container_;
};

} // namespace Opm

#endif
