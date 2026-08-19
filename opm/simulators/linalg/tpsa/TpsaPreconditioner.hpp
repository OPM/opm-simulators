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
#ifndef OPM_TPSA_PRECONDITIONER_HPP
#define OPM_TPSA_PRECONDITIONER_HPP

#include <dune/istl/operators.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/tpsa/TpsaTypes.hpp>

#include <functional>
#include <memory>
#include <type_traits>

namespace Opm
{

//! \brief Sequential operator for one of the three scalar displacement blocks.
template <typename Scalar>
using SeqDispDispOperatorT = Dune::MatrixAdapter<Linear::DispDispMatrix00T<Scalar>,
                                                 Linear::DispVector0T<Scalar>,
                                                 Linear::DispVector0T<Scalar> >;
//! \brief Sequential operator for the rotation block.
template <typename Scalar>
using SeqRotRotOperatorT = Dune::MatrixAdapter<Linear::RotRotMatrixT<Scalar>,
                                               Linear::RotVectorT<Scalar>,
                                               Linear::RotVectorT<Scalar> >;
//! \brief Sequential operator for the solid pressure block.
template <typename Scalar>
using SeqSPresSPresOperatorT = Dune::MatrixAdapter<Linear::SPresSPresMatrixT<Scalar>,
                                                   Linear::SPresVectorT<Scalar>,
                                                   Linear::SPresVectorT<Scalar> >;

#if HAVE_MPI
//! \brief Communication type of a single TPSA field.
using TpsaParComm = Dune::OwnerOverlapCopyCommunication<int, int>;

//! \brief Parallel operator for one of the three scalar displacement blocks.
template <typename Scalar>
using ParDispDispOperatorT = Dune::OverlappingSchwarzOperator<Linear::DispDispMatrix00T<Scalar>,
                                                              Linear::DispVector0T<Scalar>,
                                                              Linear::DispVector0T<Scalar>,
                                                              TpsaParComm>;
//! \brief Parallel operator for the rotation block.
template <typename Scalar>
using ParRotRotOperatorT = Dune::OverlappingSchwarzOperator<Linear::RotRotMatrixT<Scalar>,
                                                            Linear::RotVectorT<Scalar>,
                                                            Linear::RotVectorT<Scalar>,
                                                            TpsaParComm>;
//! \brief Parallel operator for the solid pressure block.
template <typename Scalar>
using ParSPresSPresOperatorT = Dune::OverlappingSchwarzOperator<Linear::SPresSPresMatrixT<Scalar>,
                                                                Linear::SPresVectorT<Scalar>,
                                                                Linear::SPresVectorT<Scalar>,
                                                                TpsaParComm>;
#endif

/*!
 * \brief Block lower-triangular preconditioner for the field-split TPSA system.
 *
 * Solves the three scalar displacement blocks, carries their contribution over
 * to the rotation and solid pressure defects, and solves those in turn.  Each
 * diagonal block gets its own Dune::FlexibleSolver, configured through the
 * `disp_disp_solver`, `rot_rot_solver` and `spres_spres_solver` sub-trees; the
 * three displacement blocks share the `disp_disp_solver` configuration.
 *
 * \tparam Scalar Field type of the matrix entries
 * \tparam DispOp Linear operator type of a single scalar displacement block,
 *                SeqDispDispOperatorT or ParDispDispOperatorT
 * \tparam RotOp Linear operator type of the rotation block, SeqRotRotOperatorT
 *               or ParRotRotOperatorT
 * \tparam SPresOp Linear operator type of the solid pressure block,
 *                 SeqSPresSPresOperatorT or ParSPresSPresOperatorT
 * \tparam Comm Communication type of the sub-solvers. Defaults to
 *              Dune::Amg::SequentialInformation, which selects the sequential
 *              implementation
 */
template <class Scalar, class DispOp, class RotOp, class SPresOp,
          class Comm = Dune::Amg::SequentialInformation>
class TpsaPreconditioner
    : public Dune::PreconditionerWithUpdate<Linear::TpsaMultiVector<Scalar>,
                                            Linear::TpsaMultiVector<Scalar> >
{
    using MultiVector = Linear::TpsaMultiVector<Scalar>;
    using DispSolver = Dune::FlexibleSolver<DispOp>;
    using RotSolver = Dune::FlexibleSolver<RotOp>;
    using SPresSolver = Dune::FlexibleSolver<SPresOp>;

public:
    //! \brief Whether this instance runs the parallel code path.
    static constexpr bool isParallel = !std::is_same_v<Comm, Dune::Amg::SequentialInformation>;

    //! \brief Field indices into the multi-vector: u_x, u_y, u_z, rotation, solid pressure.
    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;
    static constexpr auto _2 = Dune::Indices::_2;
    static constexpr auto _3 = Dune::Indices::_3;
    static constexpr auto _4 = Dune::Indices::_4;

    //! \brief The sub-solvers have no pressure equation to single out.
    static constexpr std::size_t pressureIdx = 0;

    /*!
     * \brief Sequential constructor.
     *
     * \param S View of the TPSA system matrix.
     * \param prm Configuration of the diagonal block solvers. Must hold the
     *            `disp_disp_solver`, `rot_rot_solver` and `spres_spres_solver`
     *            sub-trees
     */
    TpsaPreconditioner(const Linear::TpsaMatrixView<Scalar>& S, const PropertyTree& prm)
        requires (!isParallel);

    /*!
     * \brief Parallel constructor.
     *
     * \param S View of the TPSA system matrix.
     * \param prm Configuration of the diagonal block solvers. Must hold the
     *            `disp_disp_solver`, `rot_rot_solver` and `spres_spres_solver`
     *            sub-trees
     * \param comm Communication of a single field. All five fields live on the
     *             same dof partition, so the same object is used for each of
     *             them. Only referenced, so it must outlive the preconditioner
     */
    TpsaPreconditioner(const Linear::TpsaMatrixView<Scalar>& S,
                       const PropertyTree& prm,
                       const Comm& comm)
        requires (isParallel);

    //! \brief Nothing to prepare, the sub-solvers set themselves up.
    void pre(MultiVector&, MultiVector&) override;

    //! \brief Nothing to clean up after the last apply().
    void post(MultiVector&) override;

    //! \brief Solver category, overlapping in parallel and sequential otherwise.
    Dune::SolverCategory::Category category() const override;

    //! \brief Recompute the preconditioners of all five diagonal block solvers.
    void update() override;

    //! \brief The block solvers are rebuilt from the current matrix on update().
    bool hasPerfectUpdate() const override;

    /*!
     * \brief Apply the preconditioner, i.e. approximately solve S v = d.
     *
     * Sweeps over the fields in the order u_x, u_y, u_z, rotation, solid
     * pressure. The displacement updates are substituted into the rotation and
     * solid pressure defects before those blocks are solved.
     *
     * \param v Update, overwritten by the result
     * \param d Defect to precondition
     */
    void apply(MultiVector& v, const MultiVector& d) override;

private:
    /*!
     * \brief Build an operator and a Dune::FlexibleSolver for each diagonal block.
     *
     * \param prm Configuration of the diagonal block solvers
     */
    void initSubSolvers_(const PropertyTree& prm);

    //! \brief The system to precondition, owned by the caller.
    const Linear::TpsaMatrixView<Scalar>& S_;
    //! \brief Communication of a single field, only set in the parallel case.
    const Comm* comm_{nullptr};

    // Operators of the diagonal blocks, referenced by the solvers below.
    std::unique_ptr<DispOp> dispOp0_;
    std::unique_ptr<DispOp> dispOp1_;
    std::unique_ptr<DispOp> dispOp2_;
    std::unique_ptr<RotOp> rotOp_;
    std::unique_ptr<SPresOp> sPresOp_;

    // Solvers of the diagonal blocks.
    std::unique_ptr<DispSolver> dispSolver0_;
    std::unique_ptr<DispSolver> dispSolver1_;
    std::unique_ptr<DispSolver> dispSolver2_;
    std::unique_ptr<RotSolver> rotSolver_;
    std::unique_ptr<SPresSolver> sPresSolver_;
};

} // namespace Opm

#endif // OPM_TPSA_PRECONDITIONER_HPP
