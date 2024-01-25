/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUDILU_HPP
#define OPM_CUDILU_HPP

#include <dune/istl/preconditioner.hh>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/simulators/linalg/GraphColoring.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuMatrixDescription.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuSparseResource.hpp>
#include <vector>



namespace Opm::cuistl
{
//! \brief DILU preconditioner on the GPU.
//!
//! \tparam M The matrix type to operate on
//! \tparam X Type of the update
//! \tparam Y Type of the defect
//! \tparam l Ignored. Just there to have the same number of template arguments
//!    as other preconditioners.
//!
//! \note We assume X and Y are both CuVector<real_type>, but we leave them as template
//! arguments in case of future additions.
template <class M, class X, class Y, int l = 1>
class CuDILU : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::remove_const<M>::type;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    //! \brief Constructor.
    //!
    //!  Constructor gets all parameters to operate the prec.
    //! \param A The matrix to operate on.
    //! \param w The relaxation factor.
    //!
    explicit CuDILU(const M& A);

    //! \brief Prepare the preconditioner.
    //! \note Does nothing at the time being.
    void pre(X& x, Y& b) override;

    //! \brief Apply the preconditoner.
    void apply(X& v, const Y& d) override;

    //! \brief Post processing
    //! \note Does nothing at the moment
    void post(X& x) override;

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override;

    //! \brief Updates the matrix data.
    void update() final;


    //! \returns false
    static constexpr bool shouldCallPre()
    {
        return false;
    }

    //! \returns false
    static constexpr bool shouldCallPost()
    {
        return false;
    }


private:
    //! \brief Reference to the underlying matrix
    const M& m_cpuMatrix;
    //! \brief size_t describing the dimensions of the square block elements
    static constexpr const size_t blocksize_ = matrix_type::block_type::cols;
    //! \brief SparseTable storing each row by level
    Opm::SparseTable<size_t> m_levelSets;
    //! \brief converts from index in reordered structure to index natural ordered structure
    std::vector<int> m_reorderedToNatural;
    //! \brief converts from index in natural ordered structure to index reordered strucutre
    std::vector<int> m_naturalToReordered;
    //! \brief The A matrix stored on the gpu, and its reordred version
    CuSparseMatrix<field_type> m_gpuMatrix;
    CuSparseMatrix<field_type> m_gpuMatrixReordered;
    //! row conversion from natural to reordered matrix indices stored on the GPU
    CuVector<int> m_gpuNaturalToReorder;
    //! row conversion from reordered to natural matrix indices stored on the GPU
    CuVector<int> m_gpuReorderToNatural;
    //! \brief Stores the inverted diagonal that we use in DILU
    CuVector<field_type> m_gpuDInv;
};
} // end namespace Opm::cuistl

#endif
