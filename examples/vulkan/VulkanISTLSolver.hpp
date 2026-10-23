// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "VulkanSolver.hpp"
#include <limits>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/linalg/AbstractISTLSolver.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <sstream>

namespace Opm::Parameters
{
struct VulkanDevice {
    static constexpr auto value = "";
};
} // namespace Opm::Parameters
namespace Reslab
{
template <class TypeTag>
class VulkanISTLSolver : public Opm::AbstractISTLSolver<
                             Opm::GetPropType<TypeTag, Opm::Properties::SparseMatrixAdapter>,
                             Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>>
{
    using SparseMatrixAdapter = Opm::GetPropType<TypeTag, Opm::Properties::SparseMatrixAdapter>;
    using Vector = Opm::GetPropType<TypeTag, Opm::Properties::GlobalEqVector>;
    using Simulator = Opm::GetPropType<TypeTag, Opm::Properties::Simulator>;
    using Base = Opm::AbstractISTLSolver<SparseMatrixAdapter, Vector>;
    using Matrix = typename SparseMatrixAdapter::IstlMatrix;
    using Communication = typename Base::CommunicationType;
    using GridView = Opm::GetPropType<TypeTag, Opm::Properties::GridView>;
    using ElementContext = Opm::GetPropType<TypeTag, Opm::Properties::ElementContext>;
    using ElementChunks = Opm::ElementChunks<GridView, Dune::Partitions::All>;
    const Simulator* simulator_ = nullptr;
    std::unique_ptr<ElementChunks> chunks_;
    Opm::FlowLinearSolverParameters parameters_;
    std::unique_ptr<VulkanSolver> solver_;
    std::unique_ptr<Communication> communication_;
    Vector* rhs_ = nullptr;
    int iterations_ = 0, solveCount_ = 0;
    double solveSeconds_ = 0;
    void initialize(const Simulator& simulator)
    {
        simulator_ = &simulator;
        using namespace Opm;
        if (simulator.gridView().comm().size() != 1)
            throw std::invalid_argument("flow_vulkan currently requires one MPI process");
        if (Parameters::Get<Parameters::NonlinearSolver>() != "newton"
            || parameters_.is_nldd_local_solver_)
            throw std::invalid_argument("flow_vulkan currently supports only Newton, not NLDD");
        if (!Parameters::Get<Parameters::MatrixAddWellContributions>())
            throw std::invalid_argument(
                "flow_vulkan requires --matrix-add-well-contributions=true");
        if (parameters_.linsolver_ != "vulkan-block-jacobi"
            && parameters_.linsolver_ != "vulkan-dilu"
            && parameters_.linsolver_ != "vulkan-dilu-sweeps"
            && parameters_.linsolver_ != "vulkan-cpr")
            throw std::invalid_argument("flow_vulkan requires vulkan-block-jacobi, vulkan-dilu, "
                                        "vulkan-dilu-sweeps or vulkan-cpr");
        if (Parameters::IsSet<Parameters::LinearSolverAccelerator>()
            || parameters_.accelerator_mode_ != "none")
            throw std::invalid_argument(
                "flow_vulkan selects Vulkan through its executable; omit other accelerator flags");
        communication_ = std::make_unique<Communication>(simulator.vanguard().grid().comm());
        if (parameters_.linsolver_ == "vulkan-cpr") {
            static_assert(Opm::GetPropType<TypeTag, Opm::Properties::Indices>::pressureSwitchIdx
                          == 1);
            chunks_ = std::make_unique<ElementChunks>(
                simulator.vanguard().gridView(), Dune::Partitions::all, 1);
        }
        solver_ = std::make_unique<VulkanSolver>(Parameters::Get<Parameters::VulkanDevice>(),
                                                 parameters_.linsolver_ != "vulkan-block-jacobi",
                                                 parameters_.linsolver_ == "vulkan-dilu-sweeps"
                                                     || parameters_.linsolver_ == "vulkan-cpr",
                                                 parameters_.linsolver_ == "vulkan-cpr");
        OpmLog::info("[Vulkan] device=" + solver_->deviceName()
                     + "; FP64=true; CPU fallback=false; preconditioner=" + parameters_.linsolver_);
    }

public:
    static void registerParameters()
    {
        Opm::FlowLinearSolverParameters::registerParameters();
        Opm::Parameters::Register<Opm::Parameters::VulkanDevice>(
            "Required Vulkan hardware device name substring (CPU devices forbidden)");
    }
    explicit VulkanISTLSolver(const Simulator& simulator)
    {
        parameters_.init(simulator.vanguard().eclState().getSimulationConfig().useCPR());
        initialize(simulator);
    }
    VulkanISTLSolver(const Simulator& simulator,
                     const Opm::FlowLinearSolverParameters& parameters,
                     bool = false)
        : parameters_(parameters)
    {
        initialize(simulator);
    }
    ~VulkanISTLSolver() override
    {
        if (solver_)
            Opm::OpmLog::info("[Vulkan] completed linear solves=" + std::to_string(solveCount_)
                              + "; solve_seconds=" + std::to_string(solveSeconds_)
                              + "; CPU fallback=false");
    }
    void eraseMatrix() override
    {
        rhs_ = nullptr;
    }
    void setActiveSolver(int n) override
    {
        if (n != 0)
            throw std::invalid_argument("Only one Vulkan solver is available");
    }
    int numAvailableSolvers() const override
    {
        return 1;
    }
    void prepare(const SparseMatrixAdapter& matrix, Vector& rhs) override
    {
        prepare(matrix.istlMatrix(), rhs);
    }
    void prepare(const Matrix& matrix, Vector& rhs) override
    {
        if (matrix.N() == 0 || matrix.N() != matrix.M() || matrix.N() != rhs.size()
            || matrix[0][0].N() != 3)
            throw std::invalid_argument(
                "Vulkan requires nonempty square blocks with 3 equations per cell");
        if (matrix.N() > std::numeric_limits<uint32_t>::max() / 12
            || matrix.nonzeroes() > std::numeric_limits<uint32_t>::max() / 9)
            throw std::invalid_argument("Vulkan CSR exceeds uint32 limits");
        std::vector<uint32_t> rows;
        rows.reserve(matrix.N() * 3 + 1);
        rows.push_back(0);
        std::vector<uint32_t> columns;
        columns.reserve(matrix.nonzeroes() * 9);
        std::vector<double> values;
        values.reserve(matrix.nonzeroes() * 9);
        for (auto row = matrix.begin(); row != matrix.end(); ++row) {
            for (unsigned component = 0; component < 3; ++component) {
                for (auto entry = row->begin(); entry != row->end(); ++entry)
                    for (unsigned j = 0; j < 3; ++j) {
                        columns.push_back(entry.index() * 3 + j);
                        values.push_back((*entry)[component][j]);
                    }
                rows.push_back(columns.size());
            }
        }
        std::vector<double> pressureWeights;
        if (chunks_) {
            Vector weights(rhs.size());
            ElementContext context(*simulator_);
            Opm::Amg::getTrueImpesWeights(
                1, weights, context, simulator_->model(), *chunks_, false);
            pressureWeights.resize(3 * weights.size());
            for (size_t i = 0; i < weights.size(); ++i)
                for (unsigned j = 0; j < 3; ++j)
                    pressureWeights[3 * i + j] = weights[i][j];
        }
        solver_->prepare(std::move(rows), std::move(columns), std::move(values), pressureWeights);
        rhs_ = &rhs;
    }
    void setResidual(Vector& rhs) override
    {
        rhs_ = &rhs;
    }
    void getResidual(Vector& rhs) const override
    {
        if (!rhs_)
            throw std::logic_error("No RHS");
        rhs = *rhs_;
    }
    void setMatrix(const SparseMatrixAdapter&) override
    {
    } // prepare() owns the current matrix values.
    bool solve(Vector& x) override
    {
        if (!rhs_)
            throw std::logic_error("Vulkan solve requires prepare");
        std::vector<double> b(rhs_->size() * 3), solution;
        for (size_t i = 0; i < rhs_->size(); ++i)
            for (unsigned j = 0; j < 3; ++j)
                b[3 * i + j] = (*rhs_)[i][j];
        auto result = solver_->solve(
            b, solution, parameters_.linear_solver_reduction_, parameters_.linear_solver_maxiter_);
        iterations_ = result.iterations;
        ++solveCount_;
        solveSeconds_ += result.seconds;
        if (parameters_.linear_solver_verbosity_ > 0 || !result.converged) {
            std::ostringstream log;
            log << "[Vulkan] solve=" << solveCount_ << " iterations=" << iterations_
                << " true_relative_residual=" << result.reduction
                << " converged=" << result.converged;
            Opm::OpmLog::info(log.str());
        }
        if (!result.converged)
            throw Opm::NumericalProblem(
                "Vulkan linear solve failed true residual tolerance; no CPU fallback");
        x.resize(rhs_->size());
        for (size_t i = 0; i < x.size(); ++i)
            for (unsigned j = 0; j < 3; ++j)
                x[i][j] = solution[3 * i + j];
        return true;
    }
    int iterations() const override
    {
        return iterations_;
    }
    const Communication* comm() const override
    {
        return communication_.get();
    }
    int getSolveCount() const override
    {
        return solveCount_;
    }
};
} // namespace Reslab
