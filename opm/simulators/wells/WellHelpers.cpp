/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.
  Copyright 2020 OPM-OP AS.

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

#include <config.h>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/Well/WellInjectionControls.hpp>
#include <opm/input/eclipse/Schedule/Well/WellProductionControls.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <fmt/format.h>

#include <cstddef>
#include <vector>

namespace Opm {
namespace wellhelpers {

template<typename Scalar>
ParallelStandardWellB<Scalar>::
ParallelStandardWellB(const Matrix& B, const ParallelWellInfo& parallel_well_info)
    : B_(B), parallel_well_info_(parallel_well_info)
{}

template<typename Scalar>
template<class X, class Y>
void ParallelStandardWellB<Scalar>::
mv (const X& x, Y& y) const
{
#if !defined(NDEBUG) && HAVE_MPI
    // We need to make sure that all ranks are actually computing
    // for the same well. Doing this by checking the name of the well.
    int cstring_size = parallel_well_info_.name().size()+1;
    std::vector<int> sizes(parallel_well_info_.communication().size());
    parallel_well_info_.communication().allgather(&cstring_size, 1, sizes.data());
    std::vector<int> offsets(sizes.size()+1, 0); //last entry will be accumulated size
    std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);
    std::vector<char> cstrings(offsets[sizes.size()]);
    bool consistentWells = true;
    char* send = const_cast<char*>(parallel_well_info_.name().c_str());
    parallel_well_info_.communication().allgatherv(send, cstring_size,
                                                   cstrings.data(), sizes.data(),
                                                   offsets.data());
    for (std::size_t i = 0; i < sizes.size(); ++i)
    {
        std::string name(cstrings.data()+offsets[i]);
        if (name != parallel_well_info_.name())
        {
            if (parallel_well_info_.communication().rank() == 0)
            {
                //only one process per well logs, might not be 0 of MPI_COMM_WORLD, though
                OpmLog::error(fmt::format("Not all ranks are computing for the same well,"
                                          " should be {} but is {},", parallel_well_info_.name(), name));
            }
            consistentWells = false;
            break;
        }
    }
    parallel_well_info_.communication().barrier();
    // As not all processes are involved here we need to use MPI_Abort and hope MPI kills them all
    if (!consistentWells)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    B_.mv(x, y);

    if (this->parallel_well_info_.communication().size() > 1)
    {
        // Only do communication if we must.
        // The B matrix is basically a component-wise multiplication
        // with a vector followed by a parallel reduction. We do that
        // reduction to all ranks computing for the well to save the
        // broadcast when applying C^T.
        using YField = typename Y::block_type::value_type;
        assert(y.size() == 1);
        this->parallel_well_info_.communication().template allreduce<std::plus<YField>>(y[0].container().data(),
                                                                                        y[0].container().size());
    }
}

template<typename Scalar>
template<class X, class Y>
void ParallelStandardWellB<Scalar>::
mmv (const X& x, Y& y) const
{
    if (this->parallel_well_info_.communication().size() == 1)
    {
        // Do the same thing as before. The else branch
        // produces different rounding errors and results
        // slightly different iteration counts / well curves
        B_.mmv(x, y);
    }
    else
    {
        Y temp(y);
        mv(x, temp); // includes parallel reduction
        y -= temp;
    }
}

double computeHydrostaticCorrection(const double well_ref_depth, const double vfp_ref_depth,
                                    const double rho, const double gravity)
{
    const double dh = vfp_ref_depth - well_ref_depth;
    const double dp = rho * gravity * dh;

    return dp;
}

template<typename Scalar, typename Comm>
void sumDistributedWellEntries(Dune::DynamicMatrix<Scalar>& mat,
                               Dune::DynamicVector<Scalar>& vec,
                               const Comm& comm)
{
    // DynamicMatrix does not use one contiguous array for storing the data
    // but a DynamicVector of DynamicVectors. Hence we need to copy the data
    // to contiguous memory for MPI.
    if (comm.size() == 1)
    {
        return;
    }
    std::vector<Scalar> allEntries;
    allEntries.reserve(mat.N()*mat.M()+vec.size());
    for(const auto& row: mat)
    {
        allEntries.insert(allEntries.end(), row.begin(), row.end());
    }
    allEntries.insert(allEntries.end(), vec.begin(), vec.end());
    comm.sum(allEntries.data(), allEntries.size());
    auto pos = allEntries.begin();
    auto cols = mat.cols();
    for(auto&& row: mat)
    {
        std::copy(pos, pos + cols, &(row[0]));
        pos += cols;
    }
    assert(std::size_t(allEntries.end() - pos) == vec.size());
    std::copy(pos, allEntries.end(), &(vec[0]));
}

template <class DenseMatrix>
DenseMatrix transposeDenseDynMatrix(const DenseMatrix& M)
{
    DenseMatrix tmp{M.cols(), M.rows()};
    for (std::size_t i = 0; i < M.rows(); ++i) {
        for (std::size_t j = 0; j < M.cols(); ++j) {
            tmp[j][i] = M[i][j];
        }
    }
    return tmp;
}

bool rateControlWithZeroProdTarget(const WellProductionControls& controls,
                                   const WellProducerCMode mode)
{
    switch (mode) {
        case WellProducerCMode::ORAT:
            return controls.oil_rate == 0.0;
        case WellProducerCMode::WRAT:
            return controls.water_rate == 0.0;
        case WellProducerCMode::GRAT:
            return controls.gas_rate == 0.0;
        case WellProducerCMode::LRAT:
            return controls.liquid_rate == 0.0;
        case WellProducerCMode::CRAT:
            // Unsupported, will cause exception elsewhere, treat as nonzero target here.
            return false;
        case WellProducerCMode::RESV:
            if (controls.prediction_mode) {
                return controls.resv_rate == 0.0;
            } else {
                return controls.water_rate == 0.0 && controls.oil_rate == 0.0 && controls.gas_rate == 0.0;
            }
        default:
            return false;
    }
}


bool rateControlWithZeroInjTarget(const WellInjectionControls& controls,
                                  const WellInjectorCMode mode)
{
    switch (mode) {
        case WellInjectorCMode::RATE:
            return controls.surface_rate == 0.0;
        case WellInjectorCMode::RESV:
            return controls.reservoir_rate == 0.0;
        default:
            return false;
    }
}


template class ParallelStandardWellB<double>;

template<int Dim> using Vec = Dune::BlockVector<Dune::FieldVector<double,Dim>>;
using DynVec = Dune::BlockVector<Dune::DynamicVector<double>>;

#define INSTANCE(Dim) \
    template void ParallelStandardWellB<double>::mv<Vec<Dim>,DynVec>(const Vec<Dim>&,DynVec&) const; \
    template void ParallelStandardWellB<double>::mmv<Vec<Dim>,DynVec>(const Vec<Dim>&,DynVec&) const;

INSTANCE(1)
INSTANCE(2)
INSTANCE(3)
INSTANCE(4)
INSTANCE(5)
INSTANCE(6)

using Comm = Parallel::Communication;
template void sumDistributedWellEntries<double,Comm>(Dune::DynamicMatrix<double>& mat,
                                                     Dune::DynamicVector<double>& vec,
                                                     const Comm& comm);

using DMatrix = Dune::DynamicMatrix<double>;
template DMatrix transposeDenseDynMatrix<DMatrix>(const DMatrix&);

} // namespace wellhelpers
} // namespace Opm
