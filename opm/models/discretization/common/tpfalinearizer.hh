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
 * \copydoc Opm::FvBaseLinearizer
 */
#ifndef TPFA_LINEARIZER_HH
#define TPFA_LINEARIZER_HH

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/grid/utility/SparseTable.hpp>

#include <opm/material/common/ConditionalStorage.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/simulators/linalg/exportSystem.hpp>

// TODO: fetch via typetag of another class instead of accessing directly in this class
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>

#include <cassert>
#include <cstddef>
#include <exception>   // current_exception, rethrow_exception
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <chrono>

#include <omp.h>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/common/utility/pointerArithmetic.hpp>
#if HAVE_CUDA
#include <opm/simulators/flow/SimplifiedGpuBlackOilModel.hpp>
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl_hip/MiniMatrix.hpp>
#include <opm/simulators/linalg/gpuistl_hip/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl_hip/detail/preconditionerKernels/ILU_variants_helper_kernels.hpp>
#include <opm/simulators/linalg/gpuistl_hip/detail/gpusparse_matrix_operations.hpp>
#else
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrixWrapper.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/MiniMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/preconditionerKernels/ILU_variants_helper_kernels.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpusparse_matrix_operations.hpp>
#endif
#endif

namespace Opm::Parameters {

struct SeparateSparseSourceTerms { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

// forward declarations
template<class TypeTag>
class EcfvDiscretization;

// Moved these structs out of the class to make them visible in the GPU code.
template<class Storage = std::vector<int>>
struct FullDomain
{
    Storage cells;
};

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    FullDomain<gpuistl::GpuBuffer<int>> copy_to_gpu(FullDomain<> CPUDomain)
    {
        if (CPUDomain.cells.size() == 0) {
            OPM_THROW(std::runtime_error, "Cannot copy empty full domain to GPU.");
        }
        return FullDomain<gpuistl::GpuBuffer<int>>{
            gpuistl::GpuBuffer<int>(CPUDomain.cells)
        };
    };

    FullDomain<gpuistl::GpuView<int>> make_view(FullDomain<gpuistl::GpuBuffer<int>>& buffer)
    {
        if (buffer.cells.size() == 0) {
            OPM_THROW(std::runtime_error, "Cannot make view of empty full domain buffer.");
        }
        return FullDomain<gpuistl::GpuView<int>>{
            gpuistl::make_view(buffer.cells)
        };
    };
#endif

template <class ResidualNBInfoType,class BlockType>
struct NeighborInfoStruct
{
    unsigned int neighbor;
    ResidualNBInfoType res_nbinfo;
    BlockType* matBlockAddress;

    template <class OtherBlockType>
    NeighborInfoStruct(const NeighborInfoStruct<ResidualNBInfoType,OtherBlockType>& other)
        : neighbor(other.neighbor)
        , res_nbinfo(other.res_nbinfo)
        , matBlockAddress(nullptr)
    {
        if (other.matBlockAddress) {
            matBlockAddress = reinterpret_cast<BlockType*>(other.matBlockAddress);
        }
    }

    template <class PtrType>
    NeighborInfoStruct(unsigned int n, const ResidualNBInfoType& r, PtrType ptr)
        : neighbor(n)
        , res_nbinfo(r)
        , matBlockAddress(static_cast<BlockType*>(ptr))
    {
    }

    // Add a default constructor
    NeighborInfoStruct()
        : neighbor(0)
        , res_nbinfo()
        , matBlockAddress(nullptr)
    {
    }
};

template<class Scalar, template <class> class Storage = Opm::VectorWithDefaultAllocator>
class GpuFlowProblemVerySimple
{
public:
    GpuFlowProblemVerySimple() = default;

    GpuFlowProblemVerySimple(Storage<Scalar> alpha0,
                             Storage<Scalar> alpha1,
                             Storage<Scalar> alpha2)
        : alpha0_(alpha0)
        , alpha1_(alpha1)
        , alpha2_(alpha2)
    {}

    OPM_HOST_DEVICE Scalar getAlpha(unsigned globalIndex, unsigned boundaryFaceIndex) const
    {
        assert(boundaryFaceIndex < 3 && "SOMETHING IS WRONG WITH BOUNDARYFACEINDEX");
        if (boundaryFaceIndex == 0) {
            return alpha0_[globalIndex];
        } else if (boundaryFaceIndex == 1) {
            return alpha1_[globalIndex];
        } else {
            return alpha2_[globalIndex];
        }
    }

    OPM_HOST_DEVICE Scalar getAlpha2() const { return 1.0;}

    Storage<Scalar>& alpha0() {
        return alpha0_;
    }

    Storage<Scalar>& alpha1() {
        return alpha1_;
    }

    Storage<Scalar>& alpha2() {
        return alpha2_;
    }

private:
    Storage<Scalar> alpha0_;
    Storage<Scalar> alpha1_;
    Storage<Scalar> alpha2_;
};

namespace gpuistl {
    template<class Scalar>
    GpuFlowProblemVerySimple<Scalar, gpuistl::GpuBuffer>
    copy_to_gpu(GpuFlowProblemVerySimple<Scalar, Opm::VectorWithDefaultAllocator>& cpuProblem)
    {
        return GpuFlowProblemVerySimple<Scalar, gpuistl::GpuBuffer>(
            gpuistl::GpuBuffer<Scalar>(cpuProblem.alpha0()),
            gpuistl::GpuBuffer<Scalar>(cpuProblem.alpha1()),
            gpuistl::GpuBuffer<Scalar>(cpuProblem.alpha2())
        );
    }
    
    template<class Scalar>
    GpuFlowProblemVerySimple<Scalar, gpuistl::GpuView>
    make_view(GpuFlowProblemVerySimple<Scalar, gpuistl::GpuBuffer>& buffer)
    {
        return GpuFlowProblemVerySimple<Scalar, gpuistl::GpuView>(
            gpuistl::make_view(buffer.alpha0()),
            gpuistl::make_view(buffer.alpha1()),
            gpuistl::make_view(buffer.alpha2())
        );
    }
} // namespace gpuistl
template<class VectorBlock, class ScalarFluidState>
struct BoundaryConditionData
{
    BCType type;
    VectorBlock massRate;
    unsigned pvtRegionIdx;
    unsigned boundaryFaceIndex;
    double faceArea;
    double faceZCoord;
    ScalarFluidState exFluidState;
};

template<class BoundaryConditionData>
struct BoundaryInfo
{
    unsigned int cell;
    int dir;
    unsigned int bfIndex;
    BoundaryConditionData bcdata;
};

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
namespace  gpuistl {
    template< class MiniMatrixType, class GpuMatrixType, class CpuMatrixType, class MatrixBlockType, class ResidualNBInfoType>
    auto copy_to_gpu(const SparseTable<NeighborInfoStruct<ResidualNBInfoType, MatrixBlockType>>& cpu_neighbor_table, GpuMatrixType& gpuJacobian, CpuMatrixType& cpuJacobian)
    {
        // Convert the DUNE FieldVectors to MiniMatrix types
        using StructWithMinimatrix = NeighborInfoStruct<ResidualNBInfoType, MiniMatrixType>;
        using Scalar = typename GpuMatrixType::field_type;
        std::vector<StructWithMinimatrix> minimatrices(cpu_neighbor_table.dataSize());
        size_t idx = 0;
        for (auto e : cpu_neighbor_table.dataStorage()) {
            minimatrices[idx] = StructWithMinimatrix(e);

            Scalar* gpuBufStart = gpuJacobian.getNonZeroValues().data();
            Scalar* cpuBufStart = &(cpuJacobian[0][0][0][0]);
            Scalar* cpuPtr = &((*e.matBlockAddress)[0][0]);

            const size_t gpuNonZeroes = gpuJacobian.nonzeroes();
            const size_t cpuNonZeroes = cpuJacobian.nonzeroes();

            // To compute the length of the buffer of the cpuJacobian we here assume we have a blocked
            // BCRS matrix with square blocks and that the blocks are stored as Dune::FieldMatrix
            using CpuBlockType = typename CpuMatrixType::block_type::BaseType;

            const size_t gpuBlockSize = gpuJacobian.blockSize() * gpuJacobian.blockSize();
            const size_t cpuBlockSize = CpuBlockType::rows * CpuBlockType::cols;

            size_t gpuBufferSize = gpuNonZeroes * sizeof(typename GpuMatrixType::field_type) * gpuBlockSize;
            size_t cpuBufferSize = cpuNonZeroes * sizeof(typename CpuMatrixType::field_type) * cpuBlockSize;

            assert (gpuBufferSize == cpuBufferSize);

            // convert the pointer from CPU to GPU pointer based on offset in CPU jacobian
            Scalar* gpuPtr = ComputePtrBasedOnOffsetInOtherBuffer(
                gpuBufStart, gpuBufferSize,
                cpuBufStart, cpuBufferSize,
                cpuPtr
            );

            minimatrices[idx].matBlockAddress = reinterpret_cast<MiniMatrixType*>(gpuPtr);

            ++idx;
        }

        return SparseTable<StructWithMinimatrix, gpuistl::GpuBuffer>(
            gpuistl::GpuBuffer<StructWithMinimatrix>(minimatrices),
            gpuistl::GpuBuffer<int>(cpu_neighbor_table.rowStarts())
        );
    }

    // Handle the BoundaryInfo structs
    template<class GpuVecBlock, class GpuFluidState, class BoundaryInfoTypeGPU, class BoundaryInfoTypeCPU, typename GpuFluidSystemPtr>
    auto copy_to_gpu(const std::vector<BoundaryInfoTypeCPU>& cpu_boundary_info, GpuFluidSystemPtr* dynamicGpuFluidSystemPtr)
    {
        std::vector<BoundaryInfoTypeGPU> gpu_boundary_info;
        for (const auto& info : cpu_boundary_info) {
            gpu_boundary_info.push_back(BoundaryInfoTypeGPU{info.cell,
                                        info.dir,
                                        info.bfIndex,
                                        BoundaryConditionData<GpuVecBlock, GpuFluidState>{
                                            info.bcdata.type,
                                            GpuVecBlock(info.bcdata.massRate),
                                            info.bcdata.pvtRegionIdx,
                                            info.bcdata.boundaryFaceIndex,
                                            info.bcdata.faceArea,
                                            info.bcdata.faceZCoord,
                                            info.bcdata.exFluidState.withOtherFluidSystem(dynamicGpuFluidSystemPtr)}}); // for now we just create a dummy fluid state to avoid this weird conversion
        }

        return gpuistl::GpuBuffer<BoundaryInfoTypeGPU>(gpu_boundary_info);
    }

    // Implemented for residual_, which is a vector of FiedVectors
    // We then make a GpuBuffer of MiniVectors
    template<class CpuVecOfVecType, class GpuInnerVecType>
    auto copy_to_gpu_residual(CpuVecOfVecType& cpuVecOfVec)
    {
        std::vector<GpuInnerVecType> stdVecOfInnerVecs;
        for (const auto& vec : cpuVecOfVec) {
            stdVecOfInnerVecs.push_back((GpuInnerVecType)vec);
        }

        return gpuistl::GpuBuffer<GpuInnerVecType>(stdVecOfInnerVecs);
    }
}
#endif

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief The common code for the linearizers of non-linear systems of equations
 *
 * This class assumes that these system of equations to be linearized are stemming from
 * models that use an finite volume scheme for spatial discretization and an Euler
 * scheme for time discretization.
 */
template<class TypeTag>
class TpfaLinearizer
{
//! \cond SKIP_THIS
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;

    using Vector = GlobalEqVector;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>() };
    enum { dimWorld = GridView::dimensionworld };

    using MatrixBlock = typename SparseMatrixAdapter::MatrixBlock;
    using VectorBlock = Dune::FieldVector<Scalar, numEq>;
    using ADVectorBlock = GetPropType<TypeTag, Properties::RateVector>;

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    using MatrixBlockGPU = gpuistl::MiniMatrix<Scalar, numEq>;
    using VectorBlockGPU = gpuistl::MiniVector<Scalar, numEq>;
    using ADVectorBlockGPU = gpuistl::MiniVector<Evaluation, numEq>;
#endif

    static constexpr bool linearizeNonLocalElements =
        getPropValue<TypeTag, Properties::LinearizeNonLocalElements>();
    static constexpr bool enableFullyImplicitThermal = (getPropValue<TypeTag, Properties::EnergyModuleType>() == EnergyModules::FullyImplicitThermal);
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();
    static constexpr bool enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>();
    static const bool enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>();

    // copying the linearizer is not a good idea
    TpfaLinearizer(const TpfaLinearizer&) = delete;
//! \endcond

public:
    TpfaLinearizer()
    {
        simulatorPtr_ = nullptr;
        separateSparseSourceTerms_ = Parameters::Get<Parameters::SeparateSparseSourceTerms>();
        exportIndex_=-1;
        exportCount_=-1;
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian linearizer.
     */
    static void registerParameters()
    {
        Parameters::Register<Parameters::SeparateSparseSourceTerms>
            ("Treat well source terms all in one go, instead of on a cell by cell basis.");
    }

    /*!
     * \brief Initialize the linearizer.
     *
     * At this point we can assume that all objects in the simulator
     * have been allocated. We cannot assume that they are fully
     * initialized, though.
     *
     * \copydetails Doxygen::simulatorParam
     */
    void init(Simulator& simulator)
    {
        simulatorPtr_ = &simulator;
        eraseMatrix();
    }

    /*!
     * \brief Causes the Jacobian matrix to be recreated from scratch before the next
     *        iteration.
     *
     * This method is usally called if the sparsity pattern has changed for some
     * reason. (e.g. by modifications of the grid or changes of the auxiliary equations.)
     */
    void eraseMatrix()
    {
        jacobian_.reset();
    }

    /*!
     * \brief Linearize the full system of non-linear equations.
     *
     * The linearizationType() controls the scheme used and the focus
     * time index. The default is fully implicit scheme, and focus index
     * equal to 0, i.e. current time (end of step).
     *
     * This linearizes the spatial domain and all auxiliary equations.
     */
    void linearize()
    {
        linearizeDomain();
        linearizeAuxiliaryEquations();
    }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with the spatial domain.
     *
     * That means that the global Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by the model object.
     */
    void linearizeDomain()
    {
        int succeeded;
        try {
            linearizeDomain(fullDomain_);
            succeeded = 1;
        }
        catch (const std::exception& e) {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        catch (...) {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing"
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        succeeded = simulator_().gridView().comm().min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in linearizing the system");
        }
    }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with a part of the spatial domain.
     *
     * That means that the Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution, on the domain passed in as argument.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by the model object.
     *
     * \param domain The subdomain to linearize.
     * \param isNlddLocalSolve If true, indicates this is an NLDD local solve.
     */
    template <class SubDomainType>
    void linearizeDomain(const SubDomainType& domain, const bool isNlddLocalSolve = false)
    {
        OPM_TIMEBLOCK(linearizeDomain);
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (!jacobian_) {
            initFirstIteration_();
        }

        // Called here because it is no longer called from linearize_().
        if (isNlddLocalSolve) {
            resetSystem_(domain);
        }
        else {
            resetSystem_();
        }

        linearize_(domain, isNlddLocalSolve);
    }

    void finalize()
    { jacobian_->finalize(); }

    /*!
     * \brief Linearize the part of the non-linear system of equations that is associated
     *        with the spatial domain.
     */
    void linearizeAuxiliaryEquations()
    {
        OPM_TIMEBLOCK(linearizeAuxilaryEquations);
        // flush possible local caches into matrix structure
        jacobian_->commit();

        auto& model = model_();
        const auto& comm = simulator_().gridView().comm();
        for (unsigned auxModIdx = 0; auxModIdx < model.numAuxiliaryModules(); ++auxModIdx) {
            bool succeeded = true;
            try {
                model.auxiliaryModule(auxModIdx)->linearize(*jacobian_, residual_);
            }
            catch (const std::exception& e) {
                succeeded = false;

                std::cout << "rank " << simulator_().gridView().comm().rank()
                          << " caught an exception while linearizing:" << e.what()
                          << "\n"  << std::flush;
            }

            succeeded = comm.min(succeeded);

            if (!succeeded) {
                throw NumericalProblem("linearization of an auxiliary equation failed");
            }
        }
    }

    /*!
     * \brief Return constant reference to global Jacobian matrix backend.
     */
    const SparseMatrixAdapter& jacobian() const
    { return *jacobian_; }

    SparseMatrixAdapter& jacobian()
    { return *jacobian_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const GlobalEqVector& residual() const
    { return residual_; }

    GlobalEqVector& residual()
    { return residual_; }

    /*!
     * \brief Export block sparse linear system.
     */
    void exportSystem(const int idx, std::string& tag, const char *path="export")
    {
        const bool export_sparsity = exportIndex_ == -1;

        // increment indices and generate tag
        exportCount_ = exportIndex_ == idx ? ++exportCount_ : 0;
        exportIndex_ = idx;
        tag = fmt::format(fmt::runtime("_{:03d}_{:02d}"), exportIndex_, exportCount_);

        fmt::print(fmt::runtime("index = {:d}\n"), exportIndex_);
        fmt::print(fmt::runtime("count = {:d}\n"), exportCount_);

        Opm::exportSystem(jacobian_->istlMatrix(), residual_, export_sparsity, tag.c_str(), path);
    }

    void setLinearizationType(LinearizationType linearizationType)
    { linearizationType_ = linearizationType; }

    const LinearizationType& getLinearizationType() const
    { return linearizationType_; }

    /*!
     * \brief Return constant reference to the flowsInfo.
     *
     * (This object is only non-empty if the FLOWS keyword is true.)
     */
    const auto& getFlowsInfo() const
    { return flowsInfo_; }

    /*!
     * \brief Return constant reference to the floresInfo.
     *
     * (This object is only non-empty if the FLORES keyword is true.)
     */
    const auto& getFloresInfo() const
    { return floresInfo_; }

    /*!
     * \brief Return constant reference to the velocityInfo.
     *
     * (This object is only non-empty if the DISPERC keyword is true.)
     */
    const auto& getVelocityInfo() const
    { return velocityInfo_; }

    const auto& getNeighborInfo() const {
        return neighborInfo_;
    }


    void updateDiscretizationParameters()
    {
        updateStoredTransmissibilities();
    }

    void updateBoundaryConditionData()
    {
        for (auto& bdyInfo : boundaryInfo_) {
            const auto [type, massrateAD] = problem_().boundaryCondition(bdyInfo.cell, bdyInfo.dir);

            // Strip the unnecessary (and zero anyway) derivatives off massrate.
            VectorBlock massrate(0.0);
            for (std::size_t ii = 0; ii < massrate.size(); ++ii) {
                massrate[ii] = massrateAD[ii].value();
            }
            if (type != BCType::NONE) {
                const auto& exFluidState = problem_().boundaryFluidState(bdyInfo.cell, bdyInfo.dir);
                bdyInfo.bcdata.type = type;
                bdyInfo.bcdata.massRate = massrate;
                bdyInfo.bcdata.exFluidState = exFluidState;
            }
        }
    }

    /*!
     * \brief Returns the map of constraint degrees of freedom.
     *
     * (This object is only non-empty if the EnableConstraints property is true.)
     */
    std::map<unsigned, Constraints> constraintsMap() const
    { return {}; }

    template <class SubDomainType>
    void resetSystem_(const SubDomainType& domain)
    {
        if (!jacobian_) {
            initFirstIteration_();
        }
        for (int globI : domain.cells) {
            residual_[globI] = 0.0;
            jacobian_->clearRow(globI, 0.0);
        }
    }

private:
    Simulator& simulator_()
    { return *simulatorPtr_; }

    const Simulator& simulator_() const
    { return *simulatorPtr_; }

    Problem& problem_()
    { return simulator_().problem(); }

    const Problem& problem_() const
    { return simulator_().problem(); }

    Model& model_()
    { return simulator_().model(); }

    const Model& model_() const
    { return simulator_().model(); }

    const GridView& gridView_() const
    { return problem_().gridView(); }

    void initFirstIteration_()
    {
        // initialize the BCRS matrix for the Jacobian of the residual function
        createMatrix_();

        // initialize the Jacobian matrix and the vector for the residual function
        residual_.resize(model_().numTotalDof());
        resetSystem_();

        // initialize the sparse tables for Flows and Flores
        createFlows_();
    }

    // Construct the BCRS matrix for the Jacobian of the residual function
    void createMatrix_()
    {
        OPM_TIMEBLOCK(createMatrix);
        if (!neighborInfo_.empty()) {
            // It is ok to call this function multiple times, but it
            // should not do anything if already called.
            return;
        }
        const auto& model = model_();
        Stencil stencil(gridView_(), model_().dofMapper());

        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        using NeighborSet = std::set<unsigned>;
        std::vector<NeighborSet> sparsityPattern(model.numTotalDof());
        const Scalar gravity = problem_().gravity()[dimWorld - 1];
        unsigned numCells = model.numTotalDof();
        neighborInfo_.reserve(numCells, 6 * numCells); // Expect ~6 neighbors per cell
        std::vector<NeighborInfoCPU> loc_nbinfo;
        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                const unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                loc_nbinfo.resize(stencil.numDof() - 1); // Do not include the primary dof in neighborInfo_

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    const unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    sparsityPattern[myIdx].insert(neighborIdx);
                    if (dofIdx > 0) {
                        const Scalar trans = problem_().transmissibility(myIdx, neighborIdx);
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        const Scalar area = scvf.area();
                        const Scalar Vin = problem_().model().dofTotalVolume(myIdx);
                        const Scalar Vex = problem_().model().dofTotalVolume(neighborIdx);
                        const Scalar zIn = problem_().dofCenterDepth(myIdx);
                        const Scalar zEx = problem_().dofCenterDepth(neighborIdx);
                        const Scalar dZg = (zIn - zEx)*gravity;
                        const Scalar thpres = problem_().thresholdPressure(myIdx, neighborIdx);
                        const auto dirId = scvf.dirId();
                        auto faceDir = dirId < 0 ? FaceDir::DirEnum::Unknown
                                                 : FaceDir::FromIntersectionIndex(dirId);
                        ResidualNBInfo nbinfo{trans, area, thpres, dZg, faceDir, Vin, Vex, {}, {}, {}, {}};
                        if constexpr (enableFullyImplicitThermal) {
                            nbinfo.inAlpha = problem_().thermalHalfTransmissibility(myIdx, neighborIdx);
                            nbinfo.outAlpha = problem_().thermalHalfTransmissibility(neighborIdx, myIdx);
                        }
                        if constexpr (enableDiffusion) {
                            nbinfo.diffusivity = problem_().diffusivity(myIdx, neighborIdx);
                        }
                        if constexpr (enableDispersion) {
                            nbinfo.dispersivity = problem_().dispersivity(myIdx, neighborIdx);
                        }
                        loc_nbinfo[dofIdx - 1] = NeighborInfoCPU{neighborIdx, nbinfo, nullptr};
                    }
                }
                neighborInfo_.appendRow(loc_nbinfo.begin(), loc_nbinfo.end());
                if (problem_().nonTrivialBoundaryConditions()) {
                    for (unsigned bfIndex = 0; bfIndex < stencil.numBoundaryFaces(); ++bfIndex) {
                        const auto& bf = stencil.boundaryFace(bfIndex);
                        const int dir_id = bf.dirId();
                        // not for NNCs
                        if (dir_id < 0) {
                            continue;
                        }
                        const auto [type, massrateAD] = problem_().boundaryCondition(myIdx, dir_id);
                        // Strip the unnecessary (and zero anyway) derivatives off massrate.
                        VectorBlock massrate(0.0);
                        for (std::size_t ii = 0; ii < massrate.size(); ++ii) {
                            massrate[ii] = massrateAD[ii].value();
                        }
                        const auto& exFluidState = problem_().boundaryFluidState(myIdx, dir_id);
                        BoundaryConditionDataCPU bcdata{type,
                                                     massrate,
                                                     exFluidState.pvtRegionIndex(),
                                                     bfIndex,
                                                     bf.area(),
                                                     bf.integrationPos()[dimWorld - 1],
                                                     exFluidState};
                        boundaryInfo_.push_back({myIdx, dir_id, bfIndex, bcdata});
                    }
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        const std::size_t numAuxMod = model.numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx) {
            model.auxiliaryModule(auxModIdx)->addNeighbors(sparsityPattern);
        }

        // allocate raw matrix
        jacobian_ = std::make_unique<SparseMatrixAdapter>(simulator_());
        diagMatAddress_.resize(numCells);
        // create matrix structure based on sparsity pattern
        jacobian_->reserve(sparsityPattern);
        for (unsigned globI = 0; globI < numCells; globI++) {
            const auto& nbInfos = neighborInfo_[globI];
            diagMatAddress_[globI] = jacobian_->blockAddress(globI, globI);
            for (auto& nbInfo : nbInfos) {
                nbInfo.matBlockAddress = jacobian_->blockAddress(nbInfo.neighbor, globI);
            }
        }

#if HAVE_CUDA
        gpuJacobian_.reset(new gpuistl::GpuSparseMatrixWrapper<Scalar>(gpuistl::GpuSparseMatrixWrapper<Scalar>::fromMatrix(jacobian_->istlMatrix())));
        gpuBufferDiagMatAddress_.reset(new gpuistl::GpuBuffer<Scalar*>(gpuistl::detail::getDiagPtrs(*gpuJacobian_)));
#endif

        // Create dummy full domain.
        fullDomain_.cells.resize(numCells);
        std::iota(fullDomain_.cells.begin(), fullDomain_.cells.end(), 0);
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        // zero all matrix entries
        jacobian_->clear();

#if HAVE_CUDA
        gpuJacobian_->setToZero();
#endif
    }

    // Initialize the flows, flores, and velocity sparse tables
    void createFlows_()
    {
        OPM_TIMEBLOCK(createFlows);
        // If FLOWS/FLORES is set in any RPTRST in the schedule, then we initializate the sparse tables
        // For now, do the same also if any block flows are requested (TODO: only save requested cells...)
        // If DISPERC is in the deck, we initialize the sparse table here as well.
        const bool anyFlows = simulator_().problem().eclWriter().outputModule().getFlows().anyFlows();
        const auto& blockFlows = simulator_().problem().eclWriter().outputModule().getFlows().blockFlows();
        const bool isTemp = simulator_().vanguard().eclState().getSimulationConfig().isTemp();
        const bool anyFlores = simulator_().problem().eclWriter().outputModule().getFlows().anyFlores() || isTemp;
        const bool dispersionActive = simulator_().vanguard().eclState().getSimulationConfig().rock_config().dispersion();
        if (((!(anyFlows || !blockFlows.empty()) || !flowsInfo_.empty()) && (!anyFlores || !floresInfo_.empty())) && (!dispersionActive && !enableBioeffects)) {
            return;
        }
        const auto& model = model_();
        const auto& nncOutput = simulator_().problem().eclWriter().getOutputNnc();
        Stencil stencil(gridView_(), model_().dofMapper());
        const unsigned numCells = model.numTotalDof();
        std::unordered_multimap<int, std::pair<int, int>> nncIndices;
        std::vector<FlowInfo> loc_flinfo;
        std::vector<VelocityInfo> loc_vlinfo;
        unsigned int nncId = 0;
        VectorBlock flow(0.0);

        // Create a nnc structure to use fast lookup
        for (unsigned nncIdx = 0; nncIdx < nncOutput.size(); ++nncIdx) {
            const int ci1 = nncOutput[nncIdx].cell1;
            const int ci2 = nncOutput[nncIdx].cell2;
            nncIndices.emplace(ci1, std::make_pair(ci2, nncIdx));
        }

        if (anyFlows) {
            flowsInfo_.reserve(numCells, 6 * numCells);
        }
        else if (!blockFlows.empty()) {
            flowsInfo_.reserve(numCells, 6 * blockFlows.size());
        }
        if (anyFlores) {
            floresInfo_.reserve(numCells, 6 * numCells);
        }
        if (dispersionActive || enableBioeffects) {
            velocityInfo_.reserve(numCells, 6 * numCells);
        }

        for (const auto& elem : elements(gridView_())) {
            stencil.update(elem);
            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                const unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);
                bool blockFlowFound = false;
                if (!blockFlows.empty()) {
                    if (std::binary_search(blockFlows.begin(), blockFlows.end(),
                                           simulator_().vanguard().cartesianIndex(myIdx))) {
                        blockFlowFound = true;
                    }
                    else {
                        flowsInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.begin());
                        if (!dispersionActive && !enableBioeffects && !anyFlores) {
                            continue;
                        }
                    }
                }
                const int numFaces = stencil.numBoundaryFaces() + stencil.numInteriorFaces();
                loc_flinfo.resize(numFaces);
                loc_vlinfo.resize(stencil.numDof() - 1);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    const unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    if (dofIdx > 0) {
                        const auto scvfIdx = dofIdx - 1;
                        const auto& scvf = stencil.interiorFace(scvfIdx);
                        int faceId = scvf.dirId();
                        const int cartMyIdx = simulator_().vanguard().cartesianIndex(myIdx);
                        const int cartNeighborIdx = simulator_().vanguard().cartesianIndex(neighborIdx);
                        const auto& range = nncIndices.equal_range(cartMyIdx);
                        for (auto it = range.first; it != range.second; ++it) {
                            if (it->second.first == cartNeighborIdx){
                                // -1 gives problem since is used for the nncInput from the deck
                                faceId = -2;
                                // the index is stored to be used for writting the outputs
                                nncId = it->second.second;
                            }
                        }
                        loc_flinfo[dofIdx - 1] = FlowInfo{faceId, flow, nncId};
                        loc_vlinfo[dofIdx - 1] = VelocityInfo{flow};
                    }
                }

                for (unsigned bdfIdx = 0; bdfIdx < stencil.numBoundaryFaces(); ++bdfIdx) {
                    const auto& scvf = stencil.boundaryFace(bdfIdx);
                    const int faceId = scvf.dirId();
                    loc_flinfo[stencil.numInteriorFaces() + bdfIdx] = FlowInfo{faceId, flow, nncId};
                }

                if (anyFlows || blockFlowFound) {
                    flowsInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
                if (anyFlores) {
                    floresInfo_.appendRow(loc_flinfo.begin(), loc_flinfo.end());
                }
                if (dispersionActive || enableBioeffects) {
                    velocityInfo_.appendRow(loc_vlinfo.begin(), loc_vlinfo.end());
                }
            }
        }
    }

public:
    template<class VectorBlockType, class MatrixBlockType, class ADVectorBlockType>
    OPM_HOST_DEVICE static void setResAndJacobi(VectorBlockType& res, MatrixBlockType& bMat, const ADVectorBlockType& resid)
    {
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            res[eqIdx] = resid[eqIdx].value();
        }

        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            for (unsigned pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                // A[dofIdx][focusDofIdx][eqIdx][pvIdx] is the partial derivative of
                // the residual function 'eqIdx' for the degree of freedom 'dofIdx'
                // with regard to the focus variable 'pvIdx' of the degree of freedom
                // 'focusDofIdx'
                bMat[eqIdx][pvIdx] = resid[eqIdx].derivative(pvIdx);
            }
        }
    }

    void updateFlowsInfo()
    {
        OPM_TIMEBLOCK(updateFlows);
        const bool enableFlows = simulator_().problem().eclWriter().outputModule().getFlows().hasFlows();
        const auto& blockFlows = simulator_().problem().eclWriter().outputModule().getFlows().blockFlows();
        // We reuse the fluxes in the TEMP option
        const bool isTemp = simulator_().vanguard().eclState().getSimulationConfig().isTemp();
        const bool enableFlores = simulator_().problem().eclWriter().outputModule().getFlows().hasFlores() || isTemp;
        if (!enableFlows && !enableFlores && blockFlows.empty()) {
            return;
        }
        const unsigned int numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; ++globI) {
            OPM_TIMEBLOCK_LOCAL(linearizationForEachCell, Subsystem::Assembly);
            const auto& nbInfos = neighborInfo_[globI];
            ADVectorBlock adres(0.0);
            ADVectorBlock darcyFlux(0.0);
            const IntensiveQuantities& intQuantsIn = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
            // Flux term.
            {
                OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachCell, Subsystem::Assembly);
                short loc = 0;
                for (const auto& nbInfo : nbInfos) {
                    OPM_TIMEBLOCK_LOCAL(fluxCalculationForEachFace, Subsystem::Assembly);
                    const unsigned globJ = nbInfo.neighbor;
                    assert(globJ != globI);
                    adres = 0.0;
                    darcyFlux = 0.0;
                    const IntensiveQuantities& intQuantsEx = model_().intensiveQuantities(globJ, /*timeIdx*/ 0);
                    LocalResidual::computeFlux(adres, darcyFlux, globI, globJ, intQuantsIn,
                                               intQuantsEx, nbInfo.res_nbinfo, problem_().moduleParams());
                    adres *= nbInfo.res_nbinfo.faceArea;
                    if (!blockFlows.empty()) {
                        if (std::binary_search(blockFlows.begin(), blockFlows.end(),
                                               simulator_().vanguard().cartesianIndex(globI))) {
                            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                                flowsInfo_[globI][loc].flow[eqIdx] = adres[eqIdx].value();
                            }
                        }
                    }
                    else if (enableFlows) {
                        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                            flowsInfo_[globI][loc].flow[eqIdx] = adres[eqIdx].value();
                        }
                    }
                    if (enableFlores) {
                        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                            floresInfo_[globI][loc].flow[eqIdx] = darcyFlux[eqIdx].value();
                        }
                    }
                    ++loc;
                }
            }
        }

        // Boundary terms. Only looping over cells with nontrivial bcs.
        for (const auto& bdyInfo : boundaryInfo_) {
            if (bdyInfo.bcdata.type == BCType::NONE) {
                continue;
            }

            ADVectorBlock adres(0.0);
            const unsigned globI = bdyInfo.cell;
            const auto& nbInfos = neighborInfo_[globI];
            const IntensiveQuantities& insideIntQuants = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
            LocalResidual::computeBoundaryFlux(adres, problem_(), bdyInfo.bcdata, insideIntQuants, globI);
            adres *= bdyInfo.bcdata.faceArea;
            const unsigned bfIndex = bdyInfo.bfIndex;
            if (enableFlows) {
                for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    flowsInfo_[globI][nbInfos.size() + bfIndex].flow[eqIdx] = adres[eqIdx].value();
                }
            }
            // TODO also store Flores?
        }
    }

private:
    template <class SubDomainType>
    void linearize_(const SubDomainType& domain, bool isNlddLocalSolve)
    {
#if HAVE_CUDA
        // Make sure we have can have the domain on the GPU.
        if constexpr (std::is_same_v<SubDomainType, FullDomain<>>) {

            const bool run_assembly_on_gpu = true;

            auto enter_function = std::chrono::high_resolution_clock::now();
            /*
                One of the things I must be careful with on the GPU is all of the pointers in this class
                Block* diagMatAddress_ for instance will be pointing to CPU memory, and many helper objects
                of this nature must be converted in some way to be usable on the GPU.

                Practically I think is will work to just do pointer arithmetic and assume the same layout is used
                I think this might not be completely guaranteed though, so maybe I should be more cautious.
            */

            // This check should be removed once this is addressed by
            // for example storing the previous timesteps' values for
            // rsmax (for DRSDT) and similar.
            if (!problem_().recycleFirstIterationStorage()) {
                if (!model_().storeIntensiveQuantities() && !model_().enableStorageCache()) {
                    OPM_THROW(std::runtime_error, "Must have cached either IQs or storage when we cannot recycle.");
                }
            }

            const double dt = simulator_().timeStepSize();

            OPM_TIMEBLOCK(linearize);

            // We do not call resetSystem_() here, since that will set
            // the full system to zero, not just our part.
            // Instead, that must be called before starting the linearization.
            const bool dispersionActive = simulator_().vanguard().eclState().getSimulationConfig().rock_config().dispersion();
            const unsigned int numCells = domain.cells.size();
            const bool on_full_domain = (numCells == model_().numTotalDof());

            if (run_assembly_on_gpu) {
                hipDeviceSynchronize();
                auto prep_domain_start = std::chrono::high_resolution_clock::now();
                auto domain_buffer = copy_to_gpu(domain);
                auto domain_view = make_view(domain_buffer);
                auto prep_domain_end = std::chrono::high_resolution_clock::now();
                auto prep_domain_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_domain_end - prep_domain_start);
                std::cout << "GPU domain prep time: " << prep_domain_duration.count() << " microseconds" << std::endl;
                
                auto prep_neighbor_start = std::chrono::high_resolution_clock::now();
                auto neighborInfo_buffer = gpuistl::copy_to_gpu<MatrixBlockGPU>(neighborInfo_, *gpuJacobian_, jacobian_->istlMatrix());
                auto neighborInfo_view = gpuistl::make_view(neighborInfo_buffer);
                auto prep_neighbor_end = std::chrono::high_resolution_clock::now();
                auto prep_neighbor_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_neighbor_end - prep_neighbor_start);
                std::cout << "GPU neighborInfo prep time: " << prep_neighbor_duration.count() << " microseconds" << std::endl;

                using NeighborInfoGPU = NeighborInfoStruct<ResidualNBInfo, MatrixBlockGPU>;

                // Verify that neighborInfo_view is indeed a SparseTable with GpuView
                static_assert(std::is_same_v<decltype(neighborInfo_view),
                                            SparseTable<NeighborInfoGPU, gpuistl::GpuView>>);

                static_assert(std::is_same_v<decltype(neighborInfo_view.rowStarts()),
                                            const gpuistl::GpuView<int>&>);

                static_assert(std::is_same_v<decltype(neighborInfo_view.dataStorage()),
                                            const gpuistl::GpuView<NeighborInfoGPU>&>);

                auto diagMatAddressView = gpuistl::make_view(*gpuBufferDiagMatAddress_);

                // Take the residual_ and move it to the GPU
                // This requires going from doubles, to a blocked vector
                // This is done using a GpuBuffer which contains MiniVectors

                auto prep_residual_start = std::chrono::high_resolution_clock::now();
                auto gpuResidualBuffer = gpuistl::copy_to_gpu_residual<GlobalEqVector, VectorBlockGPU>(residual_);
                auto gpuResidualView = gpuistl::make_view(gpuResidualBuffer);
                auto prep_residual_end = std::chrono::high_resolution_clock::now();
                auto prep_residual_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_residual_end - prep_residual_start);
                std::cout << "GPU residual prep time: " << prep_residual_duration.count() << " microseconds" << std::endl;

                using CorrectTypeTagView = typename ::Opm::Properties::TTag::to_gpu_type_t<TypeTag, gpuistl::GpuView>;
                using GPUBOIQ = BlackOilIntensiveQuantities<CorrectTypeTagView>;

                using LocalResidualGPU = BlackOilLocalResidualTPFA<CorrectTypeTagView>;

                // test creating a FluidSystem that is suitable for GPU use
                auto prep_fsys_start = std::chrono::high_resolution_clock::now();
                auto& dynamicFluidSystem = FluidSystem::getNonStaticInstance();
                auto dynamicGpuFluidSystemBuffer = ::Opm::gpuistl::copy_to_gpu(dynamicFluidSystem);
                auto dynamicGpuFluidSystemView = ::Opm::gpuistl::make_view(dynamicGpuFluidSystemBuffer);
                auto prep_fsys_end = std::chrono::high_resolution_clock::now();
                auto prep_fsys_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_fsys_end - prep_fsys_start);
                std::cout << "GPU fluid system prep time: " << prep_fsys_duration.count() << " microseconds" << std::endl;

                auto prep_volumes_start = std::chrono::high_resolution_clock::now();
                std::vector<Scalar> volumes(numCells);
                for (unsigned i = 0; i < numCells; ++i) {
                    volumes[domain.cells[i]] = model_().dofTotalVolume(domain.cells[i]);
                }
                auto gpuVolumesBuffer = gpuistl::GpuBuffer<Scalar>(volumes);
                auto gpuVolumesView = gpuistl::make_view(gpuVolumesBuffer);
                auto prep_volumes_end = std::chrono::high_resolution_clock::now();
                auto prep_volumes_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_volumes_end - prep_volumes_start);
                std::cout << "GPU volumes prep time: " << prep_volumes_duration.count() << " microseconds" << std::endl;

                // We need to have a pointer to the fluidysystem that can be used inside a GPU kernel
                // Having a pointer to the view is not good enough as the view exists on the host, so
                // allocat a view on the GPU with a pointer to it via the make_gpu_shared_ptr function
                auto dynamicGpuFluidSystemPtr = gpuistl::make_gpu_shared_ptr(dynamicGpuFluidSystemView);


                using GpuScalarFluidState = typename GPUBOIQ::ScalarFluidState;
                using BoundaryConditionDataGPU = BoundaryConditionData<VectorBlockGPU, GpuScalarFluidState>;
                using BoundaryInfoGPU = BoundaryInfo<BoundaryConditionDataGPU>;
                // Copy boundary info to GPU
                gpuistl::GpuBuffer<BoundaryInfoGPU> boundaryInfo_buffer = gpuistl::copy_to_gpu<VectorBlockGPU, GpuScalarFluidState, BoundaryInfoGPU>(boundaryInfo_, dynamicGpuFluidSystemPtr.get());
                auto boundaryInfo_view = gpuistl::make_view(boundaryInfo_buffer);

                auto prep_model_start = std::chrono::high_resolution_clock::now();
                using GpuModel = GetPropType<TypeTag, Properties::GpuFIBlackOilModel>;
                GpuModel gpuModel(model_().allIntensiveQuantities0(), model_().allIntensiveQuantities1(), problem_().moduleParams());
                auto gpuModelBuffer = gpuistl::copy_to_gpu_just_find_me<TypeTag>(gpuModel, dynamicGpuFluidSystemPtr.get());
                auto gpuModelView = gpuistl::make_view_just_find_me(gpuModelBuffer);
                auto prep_model_end = std::chrono::high_resolution_clock::now();
                auto prep_model_duration = std::chrono::duration_cast<std::chrono::microseconds>(prep_model_end - prep_model_start);
                std::cout << "GPU model prep time: " << prep_model_duration.count() << " microseconds" << std::endl;


                // This is terrible, we are probably catching a very large amount of exceptions here, how to support a map on the GPU?
                // Fetch alpha values that are needed for thermal boundary condition
                std::vector<Scalar> alpha0(numCells);
                std::vector<Scalar> alpha1(numCells);
                std::vector<Scalar> alpha2(numCells);
                for (int i = 0; i < numCells; ++i) {
                    try {
                        alpha0[i] = problem_().eclTransmissibilities().thermalHalfTransBoundary(i, 0);
                    } catch (...) {
                        alpha0[i] = 0.0;
                    }
                    try {
                        alpha1[i] = problem_().eclTransmissibilities().thermalHalfTransBoundary(i, 1);
                    } catch (...) {
                        alpha1[i] = 0.0;
                    }
                    try {
                        alpha2[i] = problem_().eclTransmissibilities().thermalHalfTransBoundary(i, 2);
                    } catch (...) {
                        alpha2[i] = 0.0;
                    }
                }
                GpuFlowProblemVerySimple<Scalar> gpuFlowProblemVerySimple(alpha0, alpha1, alpha2);
                auto gpuFlowProblemVerySimpleBuffer = gpuistl::copy_to_gpu(gpuFlowProblemVerySimple);
                auto gpuFlowProblemVerySimpleView = gpuistl::make_view(gpuFlowProblemVerySimpleBuffer);
                using GpuProblem = decltype(gpuFlowProblemVerySimpleView);

                int constexpr blockSize = 256;

                hipDeviceSynchronize();
                auto start_gpu = std::chrono::high_resolution_clock::now();
                bool constexpr use_gpu = true;
                linearize_parallelization_wrapper<use_gpu, GPUBOIQ, decltype(gpuModelView), LocalResidualGPU, VectorBlockGPU, MatrixBlockGPU, ADVectorBlockGPU>(
                    numCells/*numCells*/,
                    domain_view,
                    neighborInfo_view,
                    diagMatAddressView,
                    gpuResidualView,
                    gpuModelView,
                    dt,
                    dispersionActive,
                    enableBioeffects,
                    on_full_domain,
                    gpuVolumesView);
                if (boundaryInfo_buffer.size() > 0) {
                    linearize_kernel_bc<GPUBOIQ, decltype(gpuModelView), LocalResidualGPU, VectorBlockGPU, MatrixBlockGPU, ADVectorBlockGPU><<<((boundaryInfo_buffer.size()+blockSize - 1)/blockSize), blockSize>>>(
                        diagMatAddressView,
                        gpuResidualView,
                        boundaryInfo_view,
                        gpuModelView,
                        gpuFlowProblemVerySimpleView);
                }

                hipDeviceSynchronize();
                auto end_gpu = std::chrono::high_resolution_clock::now();
                auto gpu_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_gpu - start_gpu);
                auto gpu_prep_duration = std::chrono::duration_cast<std::chrono::microseconds>(start_gpu - enter_function);
                std::cout << "GPU pre time: " << gpu_prep_duration.count() << " microseconds" << std::endl;
                std::cout << "GPU kernel time: " << gpu_duration.count() << " microseconds" << std::endl;


                auto gpu_finalize_start = std::chrono::high_resolution_clock::now();
                // Now move the gpu residual into the cpu residual
                auto cpuResidualFromGpu = gpuResidualBuffer.asStdVector();
                std::memcpy(residual_.data(), cpuResidualFromGpu.data(), numCells * numEq * sizeof(Scalar));
                {
                    auto gpuJacobianNonZeroes = gpuJacobian_->getNonZeroValues().asStdVector();
                    auto& cpuJacobian = jacobian_->istlMatrix();

                    // Copy GPU jacobian back to CPU jacobian as a single contiguous operation
                    const size_t totalMatrixSize = cpuJacobian.nonzeroes() * numEq * numEq;
                    std::memcpy(&(cpuJacobian[0][0][0][0]), gpuJacobianNonZeroes.data(), 
                               totalMatrixSize * sizeof(Scalar));
                }
                auto gpu_finalize_end = std::chrono::high_resolution_clock::now();
                auto gpu_finalize_duration = std::chrono::duration_cast<std::chrono::microseconds>(gpu_finalize_end - gpu_finalize_start);
                std::cout << "GPU post time: " << gpu_finalize_duration.count() << " microseconds" << std::endl;
            } else {
                int actual_threads = 0;

                #pragma omp parallel
                {
                    #pragma omp single
                    {
                        actual_threads = omp_get_num_threads();
                    }
                }

                printf("Using %d OpenMP threads for CPU linearization.\n", actual_threads);

                auto start_cpu = std::chrono::high_resolution_clock::now();
                linearize_parallelization_wrapper<false, IntensiveQuantities, Model, LocalResidual, VectorBlock, MatrixBlock, ADVectorBlock>(
                    numCells/*numCells*/,
                    domain,
                    neighborInfo_,
                    diagMatAddress_,
                    residual_,
                    model_(),
                    dt,
                    dispersionActive,
                    enableBioeffects,
                    on_full_domain);

                linearize_kernel_CPU_boundary<IntensiveQuantities, Model, LocalResidual, VectorBlock, MatrixBlock, ADVectorBlock>(
                    diagMatAddress_,
                    residual_,
                    boundaryInfo_);
                auto end_cpu = std::chrono::high_resolution_clock::now();
                auto cpu_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_cpu - start_cpu);

                std::cout << "CPU kernel time: " << cpu_duration.count() << " microseconds" << std::endl;
            }

            {
                // TODO: optimize the source term thing, why do we go over every cell to ask (am I a source) when we typically have 0 to a few sources?
                #pragma omp parallel for
                for (unsigned ii = 0; ii < numCells; ++ii) {

                    OPM_TIMEBLOCK_LOCAL(linearizationForEachCell, Subsystem::Assembly);
                    const unsigned globI = domain.cells[ii];
                    const auto& nbInfos = neighborInfo_[globI];

                    VectorBlock res(0.0);
                    MatrixBlock bMat(0.0);
                    ADVectorBlock adres(0.0);
                    ADVectorBlock darcyFlux(0.0);

                    const IntensiveQuantities& intQuantsIn = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
                    const Scalar volume = model_().dofTotalVolume(globI);
                    res = 0.0;
                    bMat = 0.0;
                    adres = 0.0;
                    if (separateSparseSourceTerms_) {
                        LocalResidual::computeSourceDense(adres, problem_(), intQuantsIn, globI, 0);
                    }
                    else {
                        LocalResidual::computeSource(adres, problem_(), intQuantsIn, globI, 0);
                    }
                    adres *= -volume;
                    setResAndJacobi(res, bMat, adres);
                    residual_[globI] += res;
                    //SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
                    *diagMatAddress_[globI] += bMat;
                }

                // Add sparse source terms. For now only wells.
                if (separateSparseSourceTerms_) {
                    problem_().wellModel().addReservoirSourceTerms(residual_, diagMatAddress_);
                }
            }
        }
        else {
            assert(false && "Only FullDomain is supported on GPU");
        }
#else
        printf("Wrong Linearize path HAVE_CUDA is false\n");
#endif
    }

    template<bool useGPU, class LocalIntensiveQuantities, class LocalModelClass, class LocalResidualKernel, class VectorBlockType, class MatrixBlockType, class ADVectorBlockType, class DiagPtrType, class DomainType, class NeighborSparseTable, class ResidualType>
    void linearize_parallelization_wrapper(
        const unsigned int numCells,
        const DomainType& localDomain,
        const NeighborSparseTable& localNeighborInfo,
        DiagPtrType& localDiagMatAddress,
        ResidualType& localResidual,
        LocalModelClass& localModel,
        Scalar locDT,
        bool dispersionActive,
        bool enableBioeffects,
        bool onFullDomain,
        const gpuistl::GpuView<Scalar>& localVolumes = gpuistl::GpuView<Scalar>())
    {
        if constexpr (useGPU) {
            assert(!enableBioeffects && "Bioeffects not yet supported on GPU");
            assert(!dispersionActive && "Dispersion not yet supported on GPU");
            int constexpr blockSize = 256;
            gpu_parallelize_linearization_kernel<LocalIntensiveQuantities, LocalModelClass, LocalResidualKernel, VectorBlockType, MatrixBlockType, ADVectorBlockType, DiagPtrType, DomainType, NeighborSparseTable, ResidualType><<<((numCells + blockSize - 1) / blockSize), blockSize>>>(
                numCells,
                localDomain,
                localNeighborInfo,
                localDiagMatAddress,
                localResidual,
                localModel,
                locDT,
                dispersionActive,
                enableBioeffects,
                onFullDomain,
                localVolumes);
        }
        else {
#pragma omp parallel for
            for (unsigned ii = 0; ii < numCells; ++ii) {
            linearize_kernel<false, decltype(problem_()), decltype(velocityInfo_), LocalIntensiveQuantities , LocalModelClass, LocalResidual, VectorBlockType, MatrixBlockType, ADVectorBlockType, DiagPtrType, DomainType, NeighborSparseTable, ResidualType>(
                ii,
                localDomain,
                localNeighborInfo,
                localDiagMatAddress,
                localResidual,
                localModel,
                velocityInfo_,
                locDT,
                problem_(),
                dispersionActive,
                enableBioeffects,
                onFullDomain);
            }
        }
    }

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER
    template<class LocalIntensiveQuantities, class LocalModelClass, class LocalResidualKernel, class VectorBlockType, class MatrixBlockType, class ADVectorBlockType, class DiagPtrType, class DomainType, class NeighborSparseTable, class GpuResidualView>
    __global__ static void gpu_parallelize_linearization_kernel(
        const unsigned int numCells,
        const DomainType GPU_LOCAL_domain,
        const NeighborSparseTable GPU_LOCAL_neighborInfo,
        DiagPtrType GPU_LOCAL_diagMatAddress,
        GpuResidualView GPU_LOCAL_residualView,
        LocalModelClass localModel,
        Scalar locDT,
        bool dispersionActive,
        bool enableBioeffects,
        bool onFullDomain,
        const gpuistl::GpuView<Scalar> GPU_LOCAL_volumes)
    {
        // Get the index of the cell
        const unsigned int ii = blockIdx.x * blockDim.x + threadIdx.x;

        if (ii < numCells) {
            linearize_kernel<true, int, int, LocalIntensiveQuantities, LocalModelClass, LocalResidualKernel, VectorBlockType, MatrixBlockType, ADVectorBlockType, DiagPtrType, DomainType, NeighborSparseTable, GpuResidualView>(
                ii,
                GPU_LOCAL_domain,
                GPU_LOCAL_neighborInfo,
                GPU_LOCAL_diagMatAddress,
                GPU_LOCAL_residualView,
                localModel,
                0/*dummy for velocity info*/,
                locDT,
                0/*dummy for problem type*/,
                dispersionActive,
                enableBioeffects,
                onFullDomain,
                GPU_LOCAL_volumes
                );
        }
    }
#endif

#if HAVE_CUDA && OPM_IS_COMPILING_WITH_GPU_COMPILER

    template<typename T, bool useGPU>
    using ArgType = std::conditional_t<useGPU, T, T&>;

    template<typename T, bool useGPU>
    using ConstArgType = std::conditional_t<useGPU, T, const T&>;

    template<bool useGPU, class ProblemType, class VelocityInfoType, class LocalIntensiveQuantities, class LocalModelClass, class LocalResidualKernel, class VectorBlockType, class MatrixBlockType, class ADVectorBlockType, class DiagPtrType, class DomainType, class NeighborSparseTable, class GpuResidualView>
    OPM_HOST_DEVICE static void linearize_kernel(
        const unsigned int ii,
        const ConstArgType<DomainType, useGPU> GPU_LOCAL_domain,
        const ConstArgType<NeighborSparseTable, useGPU> GPU_LOCAL_neighborInfo,
        const ConstArgType<DiagPtrType, useGPU> GPU_LOCAL_diagMatAddress,
        ArgType<GpuResidualView, useGPU> GPU_LOCAL_residualView,
        ArgType<LocalModelClass, useGPU> localModel,
        ArgType<VelocityInfoType, useGPU> localVelocityInfo,
        Scalar locDT,
        const ConstArgType<ProblemType, useGPU> localProblem,
        bool dispersionActive,
        bool enableBioeffects,
        bool on_full_domain,
        const gpuistl::GpuView<Scalar> GPU_LOCAL_volumes = gpuistl::GpuView<Scalar>()
        )
    {
        const unsigned globI = GPU_LOCAL_domain.cells[ii];
        const auto& nbInfos = GPU_LOCAL_neighborInfo[globI];
        VectorBlockType res(0.0);
        MatrixBlockType bMat(0.0);
        ADVectorBlockType adres(0.0);
        ADVectorBlockType darcyFlux(0.0);
        const LocalIntensiveQuantities& intQuantsIn = localModel.intensiveQuantities(globI, /*timeIdx*/ 0);

        // Flux term.
        {
            short loc = 0;
            for (const auto& nbInfo : nbInfos) {
                const unsigned globJ = nbInfo.neighbor;
                res = 0.0;
                bMat = 0.0;
                adres = 0.0;
                darcyFlux = 0.0;

                const LocalIntensiveQuantities& intQuantsEx = localModel.intensiveQuantities(globJ, /*timeIdx*/ 0);

                // Split because we currently do not have a gpu version of the problem object, so moduleParams are fetched from either localModel or localProblem
                if constexpr (useGPU) {
                    LocalResidualKernel::computeFlux(adres,darcyFlux, globI, globJ, intQuantsIn, intQuantsEx,
                                                    nbInfo.res_nbinfo,  localModel.moduleParams());
                } else {
                    LocalResidualKernel::computeFlux(adres,darcyFlux, globI, globJ, intQuantsIn, intQuantsEx,
                                                    nbInfo.res_nbinfo,  localProblem.moduleParams());
                }

                // currently not supported on GPU
                if constexpr (!useGPU) {
                    if (dispersionActive || enableBioeffects) {
                        for (unsigned phaseIdx = 0; phaseIdx < numEq; ++phaseIdx) {
                            localVelocityInfo[globI][loc].velocity[phaseIdx] =
                                darcyFlux[phaseIdx].value() / nbInfo.res_nbinfo.faceArea;
                        }
                        loc++;
                    }
                }

                adres *= nbInfo.res_nbinfo.faceArea;
                setResAndJacobi(res, bMat, adres);
                GPU_LOCAL_residualView[globI] += res;

                //SparseAdapter syntax:  jacobian_->addToBlock(globI, globI, bMat);
                *reinterpret_cast<MatrixBlockType*>(GPU_LOCAL_diagMatAddress[globI]) += bMat;
                bMat *= -1.0;
                //SparseAdapter syntax: jacobian_->addToBlock(globJ, globI, bMat);
                *nbInfo.matBlockAddress += bMat;
            }
        }

        // Accumulation term.
        adres = 0.0;
        {
            LocalResidualKernel::template computeStorage<Evaluation>(adres, intQuantsIn);
        }
        setResAndJacobi(res, bMat, adres);

        if constexpr (!useGPU) {
            // Either use cached storage term, or compute it on the fly.
            if (localModel.enableStorageCache()) {
                // The cached storage for timeIdx 0 (current time) is not
                // used, but after storage cache is shifted at the end of the
                // timestep, it will become cached storage for timeIdx 1.
                model_().updateCachedStorage(globI, /*timeIdx=*/0, res);
                // We should not update the storage cache here for NLDD local solves.
                // This will reset the start-of-step storage to incorrect numbers when
                // we do local solves, where the iteration number will start from 0,
                // but the starting state may not be identical to the start-of-step state.
                // Note that a full assembly must be done before local solves
                // otherwise this will be left un-updated.
                if (model_().newtonMethod().numIterations() == 0 && !isNlddLocalSolve) {
                    // Need to update the storage cache.
                    if (localProblem.recycleFirstIterationStorage()) {
                        // Assumes nothing have changed in the system which
                        // affects masses calculated from primary variables.
                        if (on_full_domain) {
                            // This is to avoid resetting the start-of-step storage
                            // to incorrect numbers when we do local solves, where the iteration
                            // number will start from 0, but the starting state may not be identical
                            // to the start-of-step state.
                            // Note that a full assembly must be done before local solves
                            // otherwise this will be left un-updated.
                            localModel.updateCachedStorage(globI, 1, res);
                        }
                    }
                    else {
                        VectorBlockType tmp;
                        const LocalIntensiveQuantities intQuantOld = localModel.intensiveQuantities(globI, 1);
                        LocalResidualKernel::template computeStorage<Scalar>(tmp, intQuantOld);
                        localModel.updateCachedStorage(globI, 1, tmp);
                    }
                }
                res -= localModel.cachedStorage(globI, 1);
            }
            else {
                OPM_TIMEBLOCK_LOCAL(computeStorage0, Subsystem::Assembly);
                VectorBlockType tmp;
                const LocalIntensiveQuantities intQuantOld = localModel.intensiveQuantities(globI, 1);
                LocalResidualKernel::template computeStorage<Scalar>(tmp, intQuantOld);
                // assume volume do not change
                res -= tmp;
            }
        } else {
            VectorBlockType tmp;
            const LocalIntensiveQuantities intQuantOld = localModel.intensiveQuantities(globI, 1);
            LocalResidualKernel::template computeStorage<Scalar>(tmp, intQuantOld);
            // assume volume do not change
            res -= tmp;
        }

        Scalar volume;
        if constexpr (useGPU) {
            volume = GPU_LOCAL_volumes[globI];
        } else {
            volume = localModel.dofTotalVolume(globI);
        }
        const Scalar storefac = volume / locDT;//volume / dt;
        res *= storefac;
        bMat *= storefac;
        GPU_LOCAL_residualView[globI] += res;
        //SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
        *reinterpret_cast<MatrixBlockType*>(GPU_LOCAL_diagMatAddress[globI]) += bMat;
    }

    template<class LocalIntensiveQuantities, class LocalModelClass, class LocalResidualKernel, class VectorBlockType, class MatrixBlockType, class ADVectorBlockType, class DiagPtrType, class GpuResidualView, class GpuBoundaryInfoView, class GpuProblem>
    __global__ static void linearize_kernel_bc(
        DiagPtrType GPU_LOCAL_diagMatAddress,
        GpuResidualView GPU_LOCAL_residualView,
        const GpuBoundaryInfoView GPU_LOCAL_boundaryInfo,
        LocalModelClass localModel,
        GpuProblem gpuProblem)
    {
        const unsigned int ii = blockIdx.x * blockDim.x + threadIdx.x;
        // Boundary terms. Only looping over cells with nontrivial bcs.
        if (ii < GPU_LOCAL_boundaryInfo.size())
        {
            // Serializing the bc handling gives correct answer but will be very slow on GPU
            // I am not sure why this did not give the correct result when using atomic adds
            // The only race condition I see is multiple threads fetching the matrix and residual
            // elements at the same time before adding something to it...
            if (GPU_LOCAL_boundaryInfo[ii].bcdata.type != BCType::NONE)
            {
                VectorBlockType res(.0);
                MatrixBlockType bMat(0.0);
                ADVectorBlockType adres(0.0);
                const unsigned globI = GPU_LOCAL_boundaryInfo[ii].cell;
                const LocalIntensiveQuantities& insideIntQuants = localModel.intensiveQuantities(globI, /*timeIdx*/ 0);
                LocalResidualKernel::computeBoundaryFlux(adres, gpuProblem, GPU_LOCAL_boundaryInfo[ii].bcdata, insideIntQuants, globI);
                adres *= GPU_LOCAL_boundaryInfo[ii].bcdata.faceArea;
                setResAndJacobi(res, bMat, adres);
                // GPU_LOCAL_residualView[globI] += res;
                auto* residualPtr = &(GPU_LOCAL_residualView.data()[globI]);
                for (int i = 0; i < numEq; ++i) {
                    atomicAdd(&((*residualPtr)[i]), res[i]);
                }
                ////SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
                // Atomic add for matrix blocks - need to add each element individually
                // *reinterpret_cast<MatrixBlockType*>(GPU_LOCAL_diagMatAddress[globI]) += bMat;
                auto* matPtr = reinterpret_cast<MatrixBlockType*>(GPU_LOCAL_diagMatAddress[globI]);
                for (int row = 0; row < bMat.size(); ++row) {
                    for (int col = 0; col < bMat.size(); ++col) {
                        Scalar* elemPtr = &((*matPtr)[row][col]);
                        atomicAdd(elemPtr, bMat[row][col]);
                    }
                }
            }
        }
    }

    template<class LocalIntensiveQuantities, class LocalModelClass, class LocalResidualKernel, class VectorBlockType, class MatrixBlockType, class ADVectorBlockType, class DiagPtrType, class GpuResidualView, class GpuBoundaryInfoView>
    void linearize_kernel_CPU_boundary(
        DiagPtrType& GPU_LOCAL_diagMatAddress,
        GpuResidualView& GPU_LOCAL_residualView,
        const GpuBoundaryInfoView& GPU_LOCAL_boundaryInfo)
    {
        // Boundary terms. Only looping over cells with nontrivial bcs.
        for (int ii = 0; ii < GPU_LOCAL_boundaryInfo.size(); ++ii)
        {
            if (GPU_LOCAL_boundaryInfo[ii].bcdata.type != BCType::NONE)
            {
                VectorBlockType res(.0);
                MatrixBlockType bMat(0.0);
                ADVectorBlockType adres(0.0);
                const unsigned globI = GPU_LOCAL_boundaryInfo[ii].cell;
                const LocalIntensiveQuantities& insideIntQuants = model_().intensiveQuantities(globI, /*timeIdx*/ 0);
                LocalResidual::computeBoundaryFlux(adres, problem_(), GPU_LOCAL_boundaryInfo[ii].bcdata, insideIntQuants, globI);
                adres *= GPU_LOCAL_boundaryInfo[ii].bcdata.faceArea;
                setResAndJacobi(res, bMat, adres);
                GPU_LOCAL_residualView[globI] += res;
                ////SparseAdapter syntax: jacobian_->addToBlock(globI, globI, bMat);
                *reinterpret_cast<MatrixBlockType*>(GPU_LOCAL_diagMatAddress[globI]) += bMat;
            }
        }
    }
#endif

    void updateStoredTransmissibilities()
    {
        if (neighborInfo_.empty()) {
            // This function was called before createMatrix_() was called.
            // We call initFirstIteration_(), not createMatrix_(), because
            // that will also initialize the residual consistently.
            initFirstIteration_();
        }

        const unsigned numCells = model_().numTotalDof();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned globI = 0; globI < numCells; globI++) {
            auto nbInfos = neighborInfo_[globI]; // nbInfos will be a SparseTable<...>::mutable_iterator_range.
            for (auto& nbInfo : nbInfos) {
                const unsigned globJ = nbInfo.neighbor;
                nbInfo.res_nbinfo.trans = problem_().transmissibility(globI, globJ);
            }
        }
    }

    Simulator* simulatorPtr_{};

    // the jacobian matrix
    std::unique_ptr<SparseMatrixAdapter> jacobian_{};
#if HAVE_CUDA
    std::unique_ptr<gpuistl::GpuSparseMatrixWrapper<Scalar>> gpuJacobian_;
    std::unique_ptr<gpuistl::GpuBuffer<Scalar*>> gpuBufferDiagMatAddress_;
#endif

    // the right-hand side
    GlobalEqVector residual_;

    LinearizationType linearizationType_{};

    using ResidualNBInfo = typename LocalResidual::ResidualNBInfo;
    using NeighborInfoCPU = NeighborInfoStruct<ResidualNBInfo, MatrixBlock>;

    SparseTable<NeighborInfoCPU> neighborInfo_{};
    std::vector<MatrixBlock*> diagMatAddress_{};

    struct FlowInfo
    {
        int faceId;
        VectorBlock flow;
        unsigned int nncId;
    };
    SparseTable<FlowInfo> flowsInfo_;
    SparseTable<FlowInfo> floresInfo_;

    struct VelocityInfo
    {
        VectorBlock velocity;
    };
    SparseTable<VelocityInfo> velocityInfo_;

    using ScalarFluidState = typename IntensiveQuantities::ScalarFluidState;

    using BoundaryConditionDataCPU = BoundaryConditionData<VectorBlock, ScalarFluidState>;

    using BoundaryInfoCPU = BoundaryInfo<BoundaryConditionDataCPU>;

    std::vector<BoundaryInfoCPU> boundaryInfo_;

    bool separateSparseSourceTerms_ = false;

    FullDomain<> fullDomain_;

    int exportIndex_;
    int exportCount_;
};
} // namespace Opm

#endif // TPFA_LINEARIZER
