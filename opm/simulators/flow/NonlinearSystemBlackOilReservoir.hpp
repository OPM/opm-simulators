/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef OPM_NONLINEAR_SYSTEM_BLACK_OIL_RESERVOIR_HEADER_INCLUDED
#define OPM_NONLINEAR_SYSTEM_BLACK_OIL_RESERVOIR_HEADER_INCLUDED

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>

#include <opm/simulators/flow/NonlinearSystem.hpp>

#include <opm/simulators/flow/BlackoilModelConvergenceMonitor.hpp>
#include <opm/simulators/flow/NonlinearSystemNldd.hpp>
#include <opm/simulators/flow/BlackoilModelProperties.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/RSTConv.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <opm/simulators/utils/ComponentName.hpp>

#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <memory>
#include <tuple>
#include <vector>

#include <fmt/format.h>

namespace Opm {

/// A model implementation for three-phase black oil.
///
/// The simulator is capable of handling three-phase problems
/// where gas can be dissolved in oil and vice versa. It
/// uses an industry-standard TPFA discretization with per-phase
/// upwind weighting of mobilities.
template <class TypeTag>
class NonlinearSystemBlackOilReservoir : public NonlinearSystem<TypeTag>
{
public:
    // ---------  Types and enums  ---------
    using ParentType = NonlinearSystem<TypeTag>;
    using Simulator = typename ParentType::Simulator;
    using Grid = typename ParentType::Grid;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalEqVector = typename ParentType::GlobalEqVector;
    using ModelParameters = BlackoilModelParameters<Scalar>;

    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

    static constexpr int numEq = Indices::numEq;
    static constexpr int contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr int contiZfracEqIdx = Indices::contiZfracEqIdx;
    static constexpr int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr int contiEnergyEqIdx = Indices::contiEnergyEqIdx;
    static constexpr int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
    static constexpr int contiFoamEqIdx = Indices::contiFoamEqIdx;
    static constexpr int contiBrineEqIdx = Indices::contiBrineEqIdx;
    static constexpr int contiMicrobialEqIdx = Indices::contiMicrobialEqIdx;
    static constexpr int contiOxygenEqIdx = Indices::contiOxygenEqIdx;
    static constexpr int contiUreaEqIdx = Indices::contiUreaEqIdx;
    static constexpr int contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
    static constexpr int contiCalciteEqIdx = Indices::contiCalciteEqIdx;
    static constexpr unsigned solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr unsigned zFractionIdx = Indices::zFractionIdx;
    static constexpr unsigned polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr unsigned polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
    static constexpr unsigned temperatureIdx = Indices::temperatureIdx;
    static constexpr unsigned foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned saltConcentrationIdx = Indices::saltConcentrationIdx;
    static constexpr unsigned microbialConcentrationIdx = Indices::microbialConcentrationIdx;
    static constexpr unsigned oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
    static constexpr unsigned ureaConcentrationIdx = Indices::ureaConcentrationIdx;
    static constexpr unsigned biofilmVolumeFractionIdx = Indices::biofilmVolumeFractionIdx;
    static constexpr unsigned calciteVolumeFractionIdx = Indices::calciteVolumeFractionIdx;

    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    using MatrixBlockType = typename SparseMatrixAdapter::MatrixBlock;
    using Mat = typename SparseMatrixAdapter::IstlMatrix;
    using BVector = Dune::BlockVector<VectorBlockType>;

    using ComponentName = ::Opm::ComponentName<FluidSystem,Indices>;

    // Helper structs
    struct CnvPvSplitData {
        std::pair<std::vector<double>, std::vector<int>> cnvPvSplit;
        std::vector<unsigned> ixCells;
    };

    struct MaxSolutionUpdateData {
        Scalar dPMax = 0.0;
        Scalar dSMax = 0.0;
        Scalar dRsMax = 0.0;
        Scalar dRvMax = 0.0;
    };

    // Output debug flags for which tolerances used.
    // NB: <windows.h> (pulled in transitively on Windows) defines STRICT as a
    // macro (=1); we also drop RELAXED defensively. This would turn "STRICT = 0"
    // into "1 = 0" in the enum below. Save the macros with push_macro and undef
    // them here; they are restored with pop_macro at the very end of this header
    // (after the _impl.hpp include), so the undef covers every in-TU use of the
    // DebugFlags::STRICT / RELAXED enumerators without leaking out to code that
    // includes this header. OPM never uses the Windows STRICT feature macro.
#if defined(_WIN32)
#  pragma push_macro("STRICT")
#  pragma push_macro("RELAXED")
#  undef STRICT
#  undef RELAXED
#endif
    enum class DebugFlags {
        STRICT = 0,
        RELAXED = 1,
        TUNINGDP = 2
    };

    // ---------  Public methods  ---------

    /// Construct the model. It will retain references to the
    /// arguments of this functions, and they are expected to
    /// remain in scope for the lifetime of the solver.
    /// \param simulator            Reference to main simulator
    /// \param[in] param            parameters
    /// \param[in] well_model       Reference to well model
    /// \param[in] terminal_output  request output to cout/cerr
    NonlinearSystemBlackOilReservoir(Simulator& simulator,
                  const ModelParameters& param,
                  BlackoilWellModel<TypeTag>& well_model,
                  const bool terminal_output);

    const EclipseState& eclState() const
    { return this->simulator_.vanguard().eclState(); }

    /// Called once before each time step.
    /// \param[in] timer                  simulation timer
    SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer);

    void initialLinearization(SimulatorReportSingle& report,
                              const int minIter,
                              const int maxIter,
                              const SimulatorTimerInterface& timer) override;

    /// Called once per nonlinear iteration.
    /// This model will perform a Newton-Raphson update, changing reservoir_state
    /// and well_state. It will also use the nonlinear_solver to do relaxation of
    /// updates if necessary.
    /// The method increases the counter in iterationContext at the end of the method.
    /// \param[in] timer                  simulation timer
    /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIteration(const SimulatorTimerInterface& timer,
                                             NonlinearSolverType& nonlinear_solver);

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationNewton(const SimulatorTimerInterface& timer,
                                                   NonlinearSolverType& nonlinear_solver);

    // compute the "relative" change of the solution between time steps
    Scalar relativeChange() const;

    /// Number of linear iterations used in last call to solveJacobianSystem().
    int linearIterationsLastSolve() const
    { return this->simulator_.model().newtonMethod().linearSolver().iterations (); }

    // Obtain reference to linear solver setup time
    double& linearSolveSetupTime()
    { return this->linear_solve_setup_time_; }

    /// Solve the Jacobian system Jx = r where J is the Jacobian and
    /// r is the residual.
    void solveJacobianSystem(BVector& x);

    /// Get solution update vector as a PrimaryVariable
    bool shouldStoreSolutionUpdate() const override;
    void prepareSolutionUpdate() override;
    void storeSolutionUpdate(const GlobalEqVector& dx) override;
    MaxSolutionUpdateData getMaxSolutionUpdate(const std::vector<unsigned>& ixCells);

    std::tuple<Scalar,Scalar>
    convergenceReduction(Parallel::Communication comm,
                         const Scalar pvSumLocal,
                         const Scalar numAquiferPvSumLocal,
                         std::vector<Scalar>& R_sum,
                         std::vector<Scalar>& maxCoeff,
                         std::vector<Scalar>& B_avg);

    /// \brief Get reservoir quantities on this process needed for convergence calculations.
    /// \return A pair of the local pore volume of interior cells and the pore volumes
    ///         of the cells associated with a numerical aquifer.
    std::pair<Scalar,Scalar>
    localConvergenceData(std::vector<Scalar>& R_sum,
                         std::vector<Scalar>& maxCoeff,
                         std::vector<Scalar>& B_avg,
                         std::vector<int>& maxCoeffCell);

    /// \brief Compute pore-volume/cell count split among "converged",
    /// "relaxed converged", "unconverged" cells based on CNV point
    /// measures. Also returns list of cells where CNV is greater than
    /// its strict tolerance
    CnvPvSplitData characteriseCnvPvSplit(const std::vector<Scalar>& B_avg, const double dt);

    /// \brief Compute the number of Newtons required by each cell in order to
    /// satisfy the solution change convergence criteria at the last time step.
    void convergencePerCell(const std::vector<Scalar>& B_avg,
                            const double dt,
                            const double tol_cnv,
                            const double tol_cnv_energy);

    ConvergenceReport
    getReservoirConvergence(const double reportTime,
                            const double dt,
                            const int maxIter,
                            std::vector<Scalar>& B_avg,
                            std::vector<Scalar>& residual_norms);

    /// Compute convergence based on total mass balance (tol_mb) and maximum
    /// residual mass balance (tol_cnv).
    /// \param[in]   timer       simulation timer
    /// \param[in]   maxIter     maximum number of iterations
    /// \param[out]  residual_norms   CNV residuals by phase
    ConvergenceReport
    getConvergence(const SimulatorTimerInterface& timer,
                   const int maxIter,
                   std::vector<Scalar>& residual_norms);

    /// Wrapper required due to not following generic API
    template<class T>
    std::vector<std::vector<Scalar> >
    computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
    { return this->computeFluidInPlace(fipnum); }

    /// Should not be called
    std::vector<std::vector<Scalar> >
    computeFluidInPlace(const std::vector<int>& /*fipnum*/) const;

    /// return the statistics of local solves accumulated for this rank
    const SimulatorReport& localAccumulatedReports() const;

    /// return the statistics of local solves accumulated for each domain on this rank
    const std::vector<SimulatorReport>& domainAccumulatedReports() const;

    /// Write the number of nonlinear iterations per cell to a file in ResInsight compatible format
    void writeNonlinearIterationsPerCell(const std::filesystem::path& odir) const;

    /// Remove the last convergence report entry and residual norms history entry.
    void popLastConvergenceReport()
    {
        this->popLastStepReport();
        this->residual_norms_history_.pop_back();
    }

    void writePartitions(const std::filesystem::path& odir) const;

    template<class FluidState, class Residual>
    void getMaxCoeff(const unsigned cell_idx,
                     const IntensiveQuantities& intQuants,
                     const FluidState& fs,
                     const Residual& modelResid,
                     const Scalar pvValue,
                     std::vector<Scalar>& B_avg,
                     std::vector<Scalar>& R_sum,
                     std::vector<Scalar>& maxCoeff,
                     std::vector<int>& maxCoeffCell);

    //! \brief Returns true if an NLDD solver exists
    bool hasNlddSolver() const
    { return nlddSolver_ != nullptr; }

protected:
    // ---------  Data members  ---------
    static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
    static constexpr bool has_extbo_ = getPropValue<TypeTag, Properties::EnableExtbo>();
    static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
    static constexpr bool has_polymermw_ = getPropValue<TypeTag, Properties::EnablePolymerMW>();
    static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal;
    static constexpr bool has_foam_ = getPropValue<TypeTag, Properties::EnableFoam>();
    static constexpr bool has_brine_ = getPropValue<TypeTag, Properties::EnableBrine>();
    static constexpr bool has_bioeffects_ = getPropValue<TypeTag, Properties::EnableBioeffects>();
    static constexpr bool has_micp_ = Indices::enableMICP;

    /// \brief The number of cells of the global grid.
    long int global_nc_;

    SolutionVector solUpd_;

    std::unique_ptr<NonlinearSystemNldd<TypeTag>> nlddSolver_; //!< Non-linear DD solver
    BlackoilModelConvergenceMonitor<Scalar> conv_monitor_;

private:
    Scalar dpMaxRel() const { return this->param_.dp_max_rel_; }
    Scalar dsMax() const { return this->param_.ds_max_; }
    Scalar drMaxRel() const { return this->param_.dr_max_rel_; }
    Scalar maxResidualAllowed() const { return this->param_.max_residual_allowed_; }
    double linear_solve_setup_time_;
    std::vector<bool> wasSwitched_;
};

} // namespace Opm

#include <opm/simulators/flow/NonlinearSystemBlackOilReservoir_impl.hpp>

#if defined(_WIN32)
#  pragma pop_macro("RELAXED")
#  pragma pop_macro("STRICT")
#endif

#endif // OPM_NONLINEAR_SYSTEM_BLACK_OIL_RESERVOIR_HEADER_INCLUDED
