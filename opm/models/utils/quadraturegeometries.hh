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
 * \copydoc Opm::QuadrialteralQuadratureGeometry
 */
#ifndef EWOMS_QUADRATURE_GEOMETRIES_HH
#define EWOMS_QUADRATURE_GEOMETRIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/type.hh>

namespace Opm {
/*!
 * \brief Quadrature geometry for quadrilaterals.
 */
template <class Scalar, unsigned dim>
class QuadrialteralQuadratureGeometry
{
public:
    enum { numCorners = (1 << dim) };

    using LocalPosition = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    Dune::GeometryType type() const
    { return Dune::GeometryType(/*topologyId=*/(1 << dim) - 1, dim); }

    template <class CornerContainer>
    void setCorners(const CornerContainer& corners, unsigned nCorners)
    {
        unsigned cornerIdx;
        for (cornerIdx = 0; cornerIdx < nCorners; ++cornerIdx) {
            for (unsigned j = 0; j < dim; ++j)
                corners_[cornerIdx][j] = corners[cornerIdx][j];
        }
        assert(cornerIdx == nCorners);

        center_ = 0;
        for (cornerIdx = 0; cornerIdx < nCorners; ++cornerIdx)
            center_ += corners_[cornerIdx];
        center_ /= nCorners;
    }

    /*!
     * \brief Returns the center of weight of the polyhedron.
     */
    const GlobalPosition& center() const
    { return center_; }

    /*!
     * \brief Convert a local coordinate into a global one.
     */
    GlobalPosition global(const LocalPosition& localPos) const
    {
        GlobalPosition globalPos(0.0);

        for (unsigned cornerIdx = 0; cornerIdx < numCorners; ++cornerIdx)
            globalPos.axpy(cornerWeight(localPos, cornerIdx),
                           corners_[cornerIdx]);

        return globalPos;
    }

    /*!
     * \brief Returns the Jacobian matrix of the local to global
     *        mapping at a given local position.
     */
    void jacobian(Dune::FieldMatrix<Scalar, dim, dim>& jac,
                  const LocalPosition& localPos) const
    {
        jac = 0.0;
        for (unsigned cornerIdx = 0; cornerIdx < numCorners; ++cornerIdx) {
            for (unsigned k = 0; k < dim; ++k) {
                Scalar dWeight_dk = (cornerIdx&  (1 << k)) ? 1 : -1;
                for (unsigned j = 0; j < dim; ++j) {
                    if (k != j) {
                        if (cornerIdx&  (1 << j))
                            dWeight_dk *= localPos[j];
                        else
                            dWeight_dk *= 1 - localPos[j];
                        ;
                    }
                }

                jac[k].axpy(dWeight_dk, corners_[cornerIdx]);
            }
        }
    }

    /*!
     * \brief Return the determinant of the Jacobian of the mapping
     *        from local to global coordinates at a given local
     *        position.
     */
    Scalar integrationElement(const LocalPosition& localPos) const
    {
        Dune::FieldMatrix<Scalar, dim, dim> jac;
        jacobian(jac, localPos);
        return jac.determinant();
    }

    /*!
     * \brief Return the position of the corner with a given index
     */
    const GlobalPosition& corner(unsigned cornerIdx) const
    { return corners_[cornerIdx]; }

    /*!
     * \brief Return the weight of an individual corner for the local
     *        to global mapping.
     */
    Scalar cornerWeight(const LocalPosition& localPos, unsigned cornerIdx) const
    {
        GlobalPosition globalPos(0.0);

        // this code is based on the Q1 finite element code from
        // dune-localfunctions
        Scalar weight = 1.0;
        for (unsigned j = 0; j < dim; ++j)
            weight *= (cornerIdx&  (1 << j)) ? localPos[j] : (1 - localPos[j]);

        return weight;
    }

private:
    GlobalPosition corners_[numCorners];
    GlobalPosition center_;
};

} // namespace Opm

#endif // EWOMS_QUADRATURE_GEOMETRY_HH
