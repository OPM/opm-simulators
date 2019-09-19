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
 * \copydoc Opm::EcfvGridCommHandleFactory
 */
#ifndef EWOMS_ECFV_GRID_COMM_HANDLE_FACTORY_HH
#define EWOMS_ECFV_GRID_COMM_HANDLE_FACTORY_HH

#include "ecfvproperties.hh"

#include <opm/models/parallel/gridcommhandles.hh>

namespace Opm {
/*!
 * \ingroup EcfvDiscretization
 *
 * \brief A class which provides types for DUNE grid handles for
 *        communication.
 *
 * This is required for parallel computations
 */
template<class TypeTag>
class EcfvGridCommHandleFactory
{
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;

public:
    /*!
     * \brief Return a handle which computes the minimum of a value
     *        for each overlapping degree of freedom across all processes.
     */
    template <class ValueType, class ArrayType>
    static std::shared_ptr<GridCommHandleMin<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> >
    minHandle(ArrayType& array, const DofMapper& dofMapper)
    {
        typedef GridCommHandleMin<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> Handle;
        return  std::shared_ptr<Handle>(new Handle(array, dofMapper));
    }

    /*!
     * \brief Return a handle which computes the maximum of a value
     *        for each overlapping degree of freedom across all processes.
     */
    template <class ValueType, class ArrayType>
    static std::shared_ptr<GridCommHandleMax<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> >
    maxHandle(ArrayType& array, const DofMapper& dofMapper)
    {
        typedef GridCommHandleMax<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> Handle;
        return  std::shared_ptr<Handle>(new Handle(array, dofMapper));
    }

    /*!
     * \brief Return a handle which computes the sum of all values
     *        all overlapping degrees of freedom across all processes.
     */
    template <class ValueType, class ArrayType>
    static std::shared_ptr<GridCommHandleSum<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> >
    sumHandle(ArrayType& array, const DofMapper& dofMapper)
    {
        typedef GridCommHandleSum<ValueType, ArrayType,  DofMapper, /*commCodim=*/0> Handle;
        return  std::shared_ptr<Handle>(new Handle(array, dofMapper));
    }
};
} // namespace Opm

#endif
