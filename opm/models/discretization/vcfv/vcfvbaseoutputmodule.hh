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
 * \copydoc Opm::VcfvBaseOutputModule
 */
#ifndef EWOMS_VCFV_VTK_BASE_OUTPUT_MODULE_HH
#define EWOMS_VCFV_VTK_BASE_OUTPUT_MODULE_HH

#include "vcfvproperties.hh"

#include <opm/models/io/baseoutputwriter.hh>

#include <string>
#include <vector>

namespace Opm {
/*!
 * \ingroup VcfvDiscretization
 *
 * \brief Implements the discretization specific parts of writing files.
 */
template<class TypeTag>
class VcfvBaseOutputModule
{
public:
    using Scalar = BaseOutputWriter::Scalar;
    using Vector = BaseOutputWriter::Vector;
    using ScalarBuffer = BaseOutputWriter::ScalarBuffer;
    using VectorBuffer = BaseOutputWriter::VectorBuffer;
    using TensorBuffer = BaseOutputWriter::TensorBuffer;

    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachScalarDofData_(BaseOutputWriter& baseWriter,
                                     ScalarBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachScalarVertexData(buffer, name.c_str()); }

    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachVectorDofData_(BaseOutputWriter& baseWriter,
                                     VectorBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachVectorVertexData(buffer, name.c_str()); }


    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    static void attachTensorDofData_(BaseOutputWriter& baseWriter,
                                     TensorBuffer& buffer,
                                     const std::string& name)
    { baseWriter.attachTensorVertexData(buffer, name.c_str()); }
};

} // namespace Opm

#endif
