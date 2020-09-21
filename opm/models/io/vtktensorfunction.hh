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
 * \copydoc Opm::VtkTensorFunction
 */
#ifndef VTK_TENSOR_FUNCTION_HH
#define VTK_TENSOR_FUNCTION_HH

#include <opm/models/io/baseoutputwriter.hh>

#include <dune/grid/io/file/vtk/function.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <string>
#include <limits>
#include <vector>

namespace Opm {

/*!
 * \brief Provides a tensor-valued function using Dune::FieldMatrix objects as elements.
 */
template <class GridView, class Mapper>
class VtkTensorFunction : public Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

    using TensorBuffer = BaseOutputWriter::TensorBuffer;

public:
    VtkTensorFunction(std::string name,
                      const GridView& gridView,
                      const Mapper& mapper,
                      const TensorBuffer& buf,
                      unsigned codim,
                      unsigned matrixColumnIdx)
        : name_(name)
        , gridView_(gridView)
        , mapper_(mapper)
        , buf_(buf)
        , codim_(codim)
        , matrixColumnIdx_(matrixColumnIdx)
    { assert(int(buf_.size()) == int(mapper_.size())); }

    virtual std::string name() const
    { return name_; }

    virtual int ncomps() const
    { return static_cast<int>(buf_[0].M()); }

    virtual double evaluate(int mycomp,
                            const Element& e,
                            const Dune::FieldVector<ctype, dim>& xi) const
    {
        size_t idx;
        if (codim_ == 0) {
            // cells. map element to the index
            idx = static_cast<size_t>(mapper_.index(e));
        }
        else if (codim_ == dim) {
            // find vertex which is closest to xi in local
            // coordinates. This code is based on Dune::P1VTKFunction
            double min = 1e100;
            int imin = -1;
            Dune::GeometryType gt = e.type();
            int n = static_cast<int>(e.subEntities(dim));
            for (int i = 0; i < n; ++i) {
                Dune::FieldVector<ctype, dim> local =
                    Dune::ReferenceElements<ctype, dim>::general(gt).position(i, dim);

                local -= xi;
                if (local.infinity_norm() < min) {
                    min = local.infinity_norm();
                    imin = i;
                }
            }

            // map vertex to an index
            idx = static_cast<size_t>(mapper_.subIndex(e, imin, codim_));
        }
        else
            throw std::logic_error("Only element and vertex based tensor fields are supported so far.");

        unsigned i = static_cast<unsigned>(mycomp);
        unsigned j = static_cast<unsigned>(matrixColumnIdx_);

        return static_cast<double>(static_cast<float>(buf_[idx][i][j]));
    }

private:
    const std::string name_;
    const GridView gridView_;
    const Mapper& mapper_;
    const TensorBuffer& buf_;
    unsigned codim_;
    unsigned matrixColumnIdx_;
};

} // namespace Opm

#endif
