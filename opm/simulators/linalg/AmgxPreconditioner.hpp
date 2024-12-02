/*
  Copyright 2024 SINTEF AS
  Copyright 2024 Equinor ASA

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

#ifndef OPM_AMGX_PRECONDITIONER_HEADER_INCLUDED
#define OPM_AMGX_PRECONDITIONER_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <amgx_c.h>

#include <vector>

namespace Amgx {

/**
 * @brief Configuration structure for AMGX parameters.
 *
 * This structure holds the configuration parameters for the AMGX solver.
 */
struct AmgxConfig {
    int determinism_flag = 0;
    int print_grid_stats = 0;
    int print_solve_stats = 0;
    std::string solver = "AMG";
    std::string algorithm = "CLASSICAL";
    std::string interpolator = "D2";
    std::string selector = "PMIS";
    std::string smoother = "BLOCK_JACOBI";
    int presweeps = 3;
    int postsweeps = 3;
    double strength_threshold = 0.5;
    int max_iters = 1;

    explicit AmgxConfig(const Opm::PropertyTree& prm) {
        determinism_flag = prm.get<int>("determinism_flag", determinism_flag);
        print_grid_stats = prm.get<int>("print_grid_stats", print_grid_stats);
        print_solve_stats = prm.get<int>("print_solve_stats", print_solve_stats);
        solver = prm.get<std::string>("solver", solver);
        algorithm = prm.get<std::string>("algorithm", algorithm);
        interpolator = prm.get<std::string>("interpolator", interpolator);
        selector = prm.get<std::string>("selector", selector);
        smoother = prm.get<std::string>("smoother", smoother);
        presweeps = prm.get<int>("presweeps", presweeps);
        postsweeps = prm.get<int>("postsweeps", postsweeps);
        strength_threshold = prm.get<double>("strength_threshold", strength_threshold);
        max_iters = prm.get<int>("max_iters", max_iters);
    }

    std::string toString() const {
        return "config_version=2, "
               "determinism_flag=" + std::to_string(determinism_flag) + ", "
               "print_grid_stats=" + std::to_string(print_grid_stats) + ", "
               "print_solve_stats=" + std::to_string(print_solve_stats) + ", "
               "solver=" + solver + ", "
               "algorithm=" + algorithm + ", "
               "interpolator=" + interpolator + ", "
               "selector=" + selector + ", "
               "smoother=" + smoother + ", "
               "presweeps=" + std::to_string(presweeps) + ", "
               "postsweeps=" + std::to_string(postsweeps) + ", "
               "strength_threshold=" + std::to_string(strength_threshold) + ", "
               "max_iters=" + std::to_string(max_iters);
    }
};

/**
 * @brief Wrapper for AMGX's AMG preconditioner.
 *
 * This class provides an interface to the AMG preconditioner from the AMGX library.
 * It is designed to work with matrices, update vectors, and defect vectors specified
 * by the template parameters.
 *
 * @tparam M The matrix type the preconditioner is for.
 * @tparam X The type of the update vector.
 * @tparam Y The type of the defect vector.
 */
template<class M, class X, class Y>
class AmgxPreconditioner : public Dune::PreconditionerWithUpdate<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for
    using matrix_type = M;
    //! \brief The domain type of the preconditioner
    using domain_type = X;
    //! \brief The range type of the preconditioner
    using range_type = Y;
    //! \brief The field type of the preconditioner
    using field_type = typename X::field_type;

    static constexpr int block_size = 1;

    /**
     * @brief Constructor for the AmgxPreconditioner class.
     *
     * Initializes the preconditioner with the given matrix and property tree.
     *
     * @param A The matrix for which the preconditioner is constructed.
     * @param prm The property tree containing configuration parameters.
     */
    AmgxPreconditioner(const M& A, const Opm::PropertyTree prm)
        : A_(A)
        , N_(A.N())
        , nnz_(A.nonzeroes())
    {
        OPM_TIMEBLOCK(prec_construct);

        // Create configuration
        AmgxConfig config(prm);
        AMGX_SAFE_CALL(AMGX_config_create(&cfg_, config.toString().c_str()));
        AMGX_SAFE_CALL(AMGX_resources_create_simple(&rsrc_, cfg_));

        // Setup frequency is set in the property tree
        setup_frequency_ = prm.get<int>("setup_frequency", 30);

        // Create solver and matrix/vector handles
        AMGX_SAFE_CALL(AMGX_solver_create(&solver_, rsrc_, AMGX_mode_dDDI, cfg_));
        AMGX_SAFE_CALL(AMGX_matrix_create(&A_amgx_, rsrc_, AMGX_mode_dDDI));
        AMGX_SAFE_CALL(AMGX_vector_create(&x_amgx_, rsrc_, AMGX_mode_dDDI));
        AMGX_SAFE_CALL(AMGX_vector_create(&b_amgx_, rsrc_, AMGX_mode_dDDI));

        // Setup matrix structure
        std::vector<int> row_ptrs(N_ + 1);
        std::vector<int> col_indices(nnz_);
        setupSparsityPattern(row_ptrs, col_indices);

        // initialize matrix with values
        const field_type* values = &(A_[0][0][0][0]);
        AMGX_SAFE_CALL(AMGX_pin_memory(const_cast<field_type*>(values), sizeof(field_type) * nnz_ * block_size * block_size));
        AMGX_SAFE_CALL(AMGX_matrix_upload_all(A_amgx_, N_, nnz_, block_size, block_size,
                                             row_ptrs.data(), col_indices.data(),
                                             values, nullptr));
        update();
    }

    /**
     * @brief Destructor for the AmgxPreconditioner class.
     *
     * Cleans up resources allocated by the preconditioner.
     */
    ~AmgxPreconditioner()
    {
        const field_type* values = &(A_[0][0][0][0]);
        AMGX_SAFE_CALL(AMGX_unpin_memory(const_cast<field_type*>(values)));
        if (solver_) {
            AMGX_SAFE_CALL(AMGX_solver_destroy(solver_));
        }
        if (x_amgx_) {
            AMGX_SAFE_CALL(AMGX_vector_destroy(x_amgx_));
        }
        if (b_amgx_) {
            AMGX_SAFE_CALL(AMGX_vector_destroy(b_amgx_));
        }
        if (A_amgx_) {
            AMGX_SAFE_CALL(AMGX_matrix_destroy(A_amgx_));
        }
        // Destroying resources and config crashes when reinitializing
        //if (rsrc_) {
        //    AMGX_SAFE_CALL(AMGX_resources_destroy(rsrc_));
        //}
        //if (cfg_) {
        //    AMGX_SAFE_CALL(AMGX_config_destroy(cfg_));
        //}
    }

    /**
     * @brief Pre-processing step before applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void pre(X& /*v*/, Y& /*d*/) override {
    }

    /**
     * @brief Applies the preconditioner to a vector.
     *
     * Performs one AMG cycle to solve the system.
     * Involves uploading vectors to AMGX, applying the preconditioner, and downloading the result.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void apply(X& v, const Y& d) override
    {
        OPM_TIMEBLOCK(prec_apply);

        // Upload vectors to AMGX
        AMGX_SAFE_CALL(AMGX_vector_upload(x_amgx_, N_, block_size, &v[0][0]));
        AMGX_SAFE_CALL(AMGX_vector_upload(b_amgx_, N_, block_size, &d[0][0]));

        // Apply preconditioner
        AMGX_SAFE_CALL(AMGX_solver_solve(solver_, b_amgx_, x_amgx_));

        // Download result
        AMGX_SAFE_CALL(AMGX_vector_download(x_amgx_, &v[0][0]));
    }

    /**
     * @brief Post-processing step after applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     */
    void post(X& /*v*/) override {
    }

    /**
     * @brief Updates the preconditioner with the current matrix values.
     *
     * This method should be called whenever the matrix values change.
     */
    void update() override
    {
        OPM_TIMEBLOCK(prec_update);
        copyMatrixToAmgx();
        if (update_counter_ == 0) {
            AMGX_SAFE_CALL(AMGX_solver_setup(solver_, A_amgx_));
        } else {
            AMGX_SAFE_CALL(AMGX_solver_resetup(solver_, A_amgx_));
        }

        ++update_counter_;
        if (update_counter_ >= setup_frequency_) {
            update_counter_ = 0;
        }
    }

    /**
     * @brief Returns the solver category.
     *
     * @return The solver category, which is sequential.
     */
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    /**
     * @brief Checks if the preconditioner has a perfect update.
     *
     * @return True, indicating that the preconditioner can be perfectly updated.
     */
    bool hasPerfectUpdate() const override
    {
        // The AMG hierarchy of the Amgx preconditioner can depend on the values of the matrix, so it must be recreated
        // when the matrix values change, at given frequency. Since this is handled internally, we return true.
        return true;
    }

private:
    /**
     * @brief Sets up the sparsity pattern for the AMGX matrix.
     *
     * This method initializes the row pointers and column indices for the AMGX matrix.
     *
     * @param row_ptrs The row pointers for the AMGX matrix.
     * @param col_indices The column indices for the AMGX matrix.
     */
    void setupSparsityPattern(std::vector<int>& row_ptrs, std::vector<int>& col_indices)
    {
        int pos = 0;
        row_ptrs[0] = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            for (auto col = row->begin(); col != row->end(); ++col) {
                col_indices[pos++] = col.index();
            }
            row_ptrs[row.index() + 1] = pos;
        }
    }

    /**
     * @brief Copies the matrix values to the AMGX matrix.
     *
     * This method updates the AMGX matrix with the current matrix values.
     * The method assumes that the sparsity structure is the same and that 
     * the values are stored in a contiguous array.
     */
    void copyMatrixToAmgx()
    {
        // Get direct pointer to matrix values
        const field_type* values = &(A_[0][0][0][0]);
        // Indexing explanation:
        // A_[0]             - First row of the matrix
        //     [0]           - First block in that row
        //        [0]        - First row within the 1x1 block
        //           [0]     - First column within the 1x1 block
        // update matrix with new values, assuming the sparsity structure is the same
        AMGX_SAFE_CALL(AMGX_matrix_replace_coefficients(A_amgx_, N_, nnz_, values, nullptr));
    }

    const M& A_; //!< The matrix for which the preconditioner is constructed.
    const int N_; //!< Number of rows in the matrix.
    const int nnz_; //!< Number of non-zero elements in the matrix.
    // Internal variables to control AMGX setup and reuse frequency
    int setup_frequency_ = -1; //!< Frequency of updating the AMG hierarchy
    int update_counter_ = 0; //!< Counter for setup updates.

    AMGX_config_handle cfg_ = nullptr; //!< The AMGX configuration handle.
    AMGX_resources_handle rsrc_ = nullptr; //!< The AMGX resources handle.
    AMGX_solver_handle solver_ = nullptr; //!< The AMGX solver handle.
    AMGX_matrix_handle A_amgx_ = nullptr; //!< The AMGX matrix handle.
    AMGX_vector_handle x_amgx_ = nullptr; //!< The AMGX solution vector handle.
    AMGX_vector_handle b_amgx_ = nullptr; //!< The AMGX right-hand side vector handle.
};

} // namespace Amgx

#endif // OPM_AMGX_PRECONDITIONER_HEADER_INCLUDED
