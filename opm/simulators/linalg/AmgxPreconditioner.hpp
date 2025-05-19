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
#include <opm/simulators/linalg/gpuistl/AmgxInterface.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <amgx_c.h>

namespace Amgx
{

using AmgxInterface = Opm::gpuistl::AmgxInterface;

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

    explicit AmgxConfig(const Opm::PropertyTree& prm)
    {
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
 * by the template parameters. The AmgxInterface class provides a unified interface
 * to AMGX's functionality, allowing for easy switching between CPU and GPU input
 * data types.
 *
 * @tparam M The matrix type
 * @tparam X The vector type for the solution
 * @tparam Y The vector type for the right-hand side
 */
template <class M, class X, class Y>
class AmgxPreconditioner : public Dune::PreconditionerWithUpdate<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for
    using matrix_type = M;
    //! \brief The field type of the matrix
    using matrix_field_type = typename M::field_type;
    //! \brief The domain type of the preconditioner
    using domain_type = X;
    //! \brief The range type of the preconditioner
    using range_type = Y;
    //! \brief The field type of the vectors
    using vector_field_type = typename X::field_type;

    static constexpr int block_size = 1;

    /**
     * @brief Constructor for AmgxPreconditioner
     *
     * Initializes the AMGX preconditioner with the given matrix and property tree.
     *
     * @param A The matrix for which to construct the preconditioner
     * @param prm The property tree containing configuration parameters
     */
    AmgxPreconditioner(const M& A, const Opm::PropertyTree prm)
        : A_(A)
        , setup_frequency_(prm.get<int>("setup_frequency", 30))
        , update_counter_(0)
    {
        OPM_TIMEBLOCK(prec_construct);

        // Create configuration
        AmgxConfig config(prm);
        cfg_ = AmgxInterface::createConfig(config.toString());
        rsrc_ = AmgxInterface::createResources(cfg_);

        // Determine appropriate AMGX mode based on matrix and vector types
        amgx_mode_ = AmgxInterface::determineAmgxMode<matrix_field_type, vector_field_type>();

        // Create solver, matrix, and vector handles, given the mode
        solver_ = AmgxInterface::createSolver(rsrc_, amgx_mode_, cfg_);
        A_amgx_ = AmgxInterface::createMatrix(rsrc_, amgx_mode_);
        x_amgx_ = AmgxInterface::createVector(rsrc_, amgx_mode_);
        b_amgx_ = AmgxInterface::createVector(rsrc_, amgx_mode_);

        // Initialize matrix structure and values
        AmgxInterface::initializeMatrix(A_, A_amgx_);

        // Initialize vectors with proper dimensions
        const int N = A_.N();
        AmgxInterface::initializeVector(N, block_size, x_amgx_);
        AmgxInterface::initializeVector(N, block_size, b_amgx_);

        // Perform initial update
        update();
    }

    /**
     * @brief Destructor for AmgxPreconditioner
     *
     * Cleans up resources allocated by the preconditioner.
     */
    ~AmgxPreconditioner()
    {
        AmgxInterface::destroySolver(solver_);
        AmgxInterface::destroyVector(x_amgx_);
        AmgxInterface::destroyVector(b_amgx_);
        AmgxInterface::destroyMatrix(A_amgx_, A_);
        AmgxInterface::destroyResources(rsrc_);
        AmgxInterface::destroyConfig(cfg_);
    }

    /**
     * @brief Pre-processing step before applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void pre(X& /*v*/, Y& /*d*/) override
    {
    }

    /**
     * @brief Applies the preconditioner to a vector.
     *
     * Performs one AMG cycle to solve the system.
     * Involves uploading vectors to AMGX, applying the preconditioner,
     * and transferring the result back to the vector.
     *
     * @param v The update vector.
     * @param d The defect vector.
     */
    void apply(X& v, const Y& d) override
    {
        OPM_TIMEBLOCK(prec_apply);

        // Transfer vectors to AMGX
        AmgxInterface::transferVectorToAmgx(v, x_amgx_);
        AmgxInterface::transferVectorToAmgx(d, b_amgx_);

        // Apply preconditioner
        AMGX_SAFE_CALL(AMGX_solver_solve(solver_, b_amgx_, x_amgx_));

        // Transfer result back
        AmgxInterface::transferVectorFromAmgx(x_amgx_, v);
    }

    /**
     * @brief Post-processing step after applying the preconditioner.
     *
     * This method is currently a no-op.
     *
     * @param v The update vector.
     */
    void post(X& /*v*/) override
    {
    }

    /**
     * @brief Updates the preconditioner with the current matrix values.
     *
     * This method should be called whenever the matrix values change.
     */
    void update() override
    {
        OPM_TIMEBLOCK(prec_update);

        AmgxInterface::updateMatrixValues(A_, A_amgx_);

        // Setup or resetup the solver based on update counter
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
     * @return The solver category (sequential).
     */
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    /**
     * @brief Checks if the preconditioner has a perfect update.
     *
     * @return True, as the preconditioner can be perfectly updated.
     */
    bool hasPerfectUpdate() const override
    {
        // The AMG hierarchy of the Amgx preconditioner can depend on the values of the matrix, so it must be recreated
        // when the matrix values change, at a given frequency. Since this is handled internally, we return true.
        return true;
    }

private:
    // Reference to the input matrix
    const M& A_;

    // AMGX handles
    AMGX_config_handle cfg_ = nullptr;
    AMGX_resources_handle rsrc_ = nullptr;
    AMGX_solver_handle solver_ = nullptr;
    AMGX_matrix_handle A_amgx_ = nullptr;
    AMGX_vector_handle x_amgx_ = nullptr;
    AMGX_vector_handle b_amgx_ = nullptr;
    AMGX_Mode amgx_mode_;

    // Frequency of updating the AMG hierarchy
    int setup_frequency_;
    // Counter for setup updates
    int update_counter_;
};

} // namespace Amgx

#endif // OPM_AMGX_PRECONDITIONER_HEADER_INCLUDED
