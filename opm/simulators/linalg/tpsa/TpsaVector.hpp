// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_TPSA_VECTOR_HPP
#define OPM_TPSA_VECTOR_HPP

#include <opm/simulators/linalg/tpsa/TpsaTypes.hpp>

#include <dune/common/fvector.hh>

#include <array>
#include <cmath>
#include <cstddef>

namespace Opm::Linear
{

/*!
 * \brief Vector type for TPSA linear elasticity
 *
 * Companion to TpsaMatrix class. This is a wrapper around Dune::MultiTypeBlockVector type to
 * distribute a 7x1 dense vector to sub-vector fields corresponding to the sub-matrix fields in
 * TpsaMatrix
 *
 * \tparam ScalarT Field type of the vector entries.
 */
template <class ScalarT>
class TpsaVector
{
public:
    //! \brief Field type of the vector entries.
    using Scalar = ScalarT;

    //! \copydoc Scalar
    using field_type = ScalarT;

    //! \brief Dense block holding the seven equations of a single cell.
    using EqVector = Dune::FieldVector<Scalar, numTpsaEq>;

    //! \copydoc EqVector
    using block_type = EqVector;

    //! \brief What the linear solver operates on.
    using IstlVector = TpsaMultiVector<Scalar>;

    //! \brief Type used for sizes and indices.
    using size_type = std::size_t;

    /*!
     * \brief Mutable reference to the seven equations of a single cell.
     *
     * Reads and writes are scattered over the five sub-vectors.
     */
    class EntryProxy
    {
    public:
        /*!
         * \brief Construct a handle on the equations of one cell.
         *
         * \param[in] v Multi-type vector holding the sub-vectors. Must outlive the proxy,
         *              which only stores a pointer to it.
         * \param[in] dofIdx Index of the cell the proxy refers to.
         */
        EntryProxy(IstlVector& v, std::size_t dofIdx)
            : v_(&v)
            , i_(dofIdx)
        {
        }

        /*!
         * \brief Number of equations the proxy exposes.
         *
         * \return Number of TPSA equations per cell.
         */
        static constexpr std::size_t size()
        {
            return numTpsaEq;
        }

        /*!
         * \brief Access one equation of the cell.
         *
         * \param[in] eqIdx Equation index in [0, numTpsaEq)
         *
         * \return Reference to the entry in the sub-vector it belongs to.
         */
        Scalar& operator[](std::size_t eqIdx)
        {
            return at_(eqIdx);
        }

        /*!
         * \brief Read one equation of the cell.
         *
         * \param[in] eqIdx Equation index in [0, numTpsaEq)
         *
         * \return Value of the entry.
         */
        Scalar operator[](std::size_t eqIdx) const
        {
            return at_(eqIdx);
        }

        /*!
         * \brief Gathers to an EqVector block vector
         *
         * \warning This function requires an explicit return type EqVector, contrary to
         * "auto b = v[i]" which returns a EntryProxy type
         *
         * \return The seven equations of the cell, gathered into a dense block.
         */
        operator EqVector() const;

        /*!
         * \brief Set all seven equations of the cell to a scalar.
         *
         * \param[in] value Value assigned to every equation.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator=(Scalar value);

        /*!
         * \brief Scatter a dense 7x1 block into the sub-vectors.
         *
         * \param[in] value Dense block the entries are copied from.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator=(const EqVector& value);

        /*!
         * \brief Copy the equations of one cell to another.
         *
         * \param[in] other Proxy the entries are copied from.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator=(const EntryProxy& other);

        /*!
         * \brief Add a dense 7x1 block to the equations of the cell.
         *
         * \param[in] value Dense block added to the stored entries.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator+=(const EqVector& value);

        /*!
         * \brief Subtract a dense 7x1 block from the equations of the cell.
         *
         * \param[in] value Dense block subtracted from the stored entries.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator-=(const EqVector& value);

        /*!
         * \brief Scale all seven equations of the cell.
         *
         * \param[in] factor Factor the stored entries are multiplied by.
         *
         * \return Reference to this proxy.
         */
        EntryProxy& operator*=(Scalar factor);

    private:
        /*!
         * \brief Map an equation index onto the entry of the sub-vector holding it.
         *
         * \param[in] eqIdx Equation index.
         *
         * \return Reference to the corresponding sub-vector entry.
         */
        Scalar& at_(std::size_t eqIdx) const;

        IstlVector* v_;
        std::size_t i_;
    };

    //! \brief Construct an empty vector; call resize() before use.
    TpsaVector() = default;

    /*!
     * \brief Construct a vector sized for the given number of cells.
     *
     * \param[in] numDof Number of degrees of freedom, i.e. cells.
     */
    explicit TpsaVector(std::size_t numDof)
    {
        resize(numDof);
    }

    /*!
     * \brief Resize every sub-vector to the given number of cells.
     *
     * \param[in] numDof Number of degrees of freedom, i.e. cells.
     */
    void resize(std::size_t numDof);

    /*!
     * \brief Number of cells the vector holds equations for.
     *
     * \return Number of degrees of freedom.
     */
    std::size_t size() const
    {
        return size_;
    }

    //! \copydoc size()
    std::size_t N() const
    {
        return size_;
    }

    /*!
     * \brief Set every entry of every field to a scalar.
     *
     * \param[in] value Value assigned to all entries.
     *
     * \return Reference to this vector.
     */
    TpsaVector& operator=(Scalar value);

    /*!
     * \brief Add another vector field by field.
     *
     * \param[in] other Vector added to this one. Must have the same size.
     *
     * \return Reference to this vector.
     */
    TpsaVector& operator+=(const TpsaVector& other);

    /*!
     * \brief Subtract another vector field by field.
     *
     * \param[in] other Vector subtracted from this one. Must have the same size.
     *
     * \return Reference to this vector.
     */
    TpsaVector& operator-=(const TpsaVector& other);

    /*!
     * \brief Scale every entry of every field.
     *
     * \param[in] factor Factor all entries are multiplied by.
     *
     * \return Reference to this vector.
     */
    TpsaVector& operator*=(Scalar factor);

    /*!
     * \brief Multiply each field by factor.
     *
     * The counterpart of TpsaMatrix::scaleFields().
     *
     * \param[in] factors Factor for each field, in field order.
     *
     * \note Scales in place.
     */
    void scaleFields(const std::array<Scalar, numTpsaFields>& factors);

    /*!
     * \brief Sum of the absolute values of all entries.
     *
     * \return The 1-norm of the vector.
     */
    Scalar one_norm() const;

    /*!
     * \brief Sum of the squares of all entries.
     *
     * \return The squared 2-norm of the vector.
     */
    Scalar two_norm2() const;

    /*!
     * \brief Euclidean norm of the vector.
     *
     * \return The 2-norm of the vector.
     */
    Scalar two_norm() const
    {
        return std::sqrt(two_norm2());
    }

    /*!
     * \brief Largest absolute value over all entries.
     *
     * \return The infinity norm of the vector.
     */
    Scalar infinity_norm() const;

    /*!
     * \brief The seven equations of cell dofIdx, gathered into a dense block.
     *
     * \param[in] dofIdx Index of the cell.
     *
     * \return Dense block holding the seven equations of that cell.
     */
    EqVector operator[](std::size_t dofIdx) const;

    /*!
     * \brief Writable handle on the seven equations of cell dofIdx.
     *
     * \param[in] dofIdx Index of the cell.
     *
     * \return Proxy scattering reads and writes over the five sub-vectors.
     */
    EntryProxy operator[](std::size_t dofIdx)
    {
        return EntryProxy(v_, dofIdx);
    }

    /*!
     * \brief The underlying multi-type vector handed to the linear solver.
     *
     * \return Reference to the multi-type block vector holding the five sub-vectors.
     */
    IstlVector& istlVector()
    {
        return v_;
    }

    //! \copydoc istlVector()
    const IstlVector& istlVector() const
    {
        return v_;
    }

private:
    IstlVector v_{};
    std::size_t size_{0};
};

} // namespace Opm::Linear

#endif // OPM_TPSA_VECTOR_HPP
