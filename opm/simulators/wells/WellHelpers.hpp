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


#ifndef OPM_WELLHELPERS_HEADER_INCLUDED
#define OPM_WELLHELPERS_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <vector>

namespace Opm {




    namespace wellhelpers
    {

        /// \brief A wrapper around the B matrix for distributed wells
        ///
        /// For standard wells the B matrix, is basically a multiplication
        /// of the equation of the perforated cells followed by a reduction
        /// (summation) of these to the well equations.
        ///
        /// This class does that in the functions mv and mmv (from the DUNE
        /// matrix interface.
        ///
        /// \tparam Scalar The scalar used for the computation.
        template<typename Scalar>
        class ParallelStandardWellB
        {
        public:
            using Block = Dune::DynamicMatrix<Scalar>;
            using Matrix = Dune::BCRSMatrix<Block>;

            ParallelStandardWellB(const Matrix& B, const ParallelWellInfo& parallel_well_info)
                : B_(&B), parallel_well_info_(&parallel_well_info)
            {}

            //! y = A x
            template<class X, class Y>
            void mv (const X& x, Y& y) const
            {
#if !defined(NDEBUG) && HAVE_MPI
                // We need to make sure that all ranks are actually computing
                // for the same well. Doing this by checking the name of the well.
                int cstring_size = parallel_well_info_->name().size()+1;
                std::vector<int> sizes(parallel_well_info_->communication().size());
                parallel_well_info_->communication().allgather(&cstring_size, 1, sizes.data());
                std::vector<int> offsets(sizes.size()+1, 0); //last entry will be accumulated size
                std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);
                std::vector<char> cstrings(offsets[sizes.size()]);
                bool consistentWells = true;
                char* send = const_cast<char*>(parallel_well_info_->name().c_str());
                parallel_well_info_->communication().allgatherv(send, cstring_size,
                                                   cstrings.data(), sizes.data(),
                                                   offsets.data());
                for(std::size_t i = 0; i < sizes.size(); ++i)
                {
                    std::string name(cstrings.data()+offsets[i]);
                    if (name != parallel_well_info_->name())
                    {
                        if (parallel_well_info_->communication().rank() == 0)
                        {
                            //only one process per well logs, might not be 0 of MPI_COMM_WORLD, though
                            std::string msg = std::string("Fatal Error: Not all ranks are computing for the same well")
                                          + " well should be " + parallel_well_info_->name() + " but is "
                                + name;
                            OpmLog::debug(msg);
                        }
                        consistentWells = false;
                        break;
                    }
                }
                parallel_well_info_->communication().barrier();
                // As not all processes are involved here we need to use MPI_Abort and hope MPI kills them all
                if (!consistentWells)
                {
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
#endif
                B_->mv(x, y);

                if (this->parallel_well_info_->communication().size() > 1)
                {
                    // Only do communication if we must.
                    // The B matrix is basically a component-wise multiplication
                    // with a vector followed by a parallel reduction. We do that
                    // reduction to all ranks computing for the well to save the
                    // broadcast when applying C^T.
                    using YField = typename Y::block_type::value_type;
                    assert(y.size() == 1);
                    this->parallel_well_info_->communication().template allreduce<std::plus<YField>>(y[0].container().data(),
                                                                                                     y[0].container().size());
                }
            }

            //! y = A x
            template<class X, class Y>
            void mmv (const X& x, Y& y) const
            {
                if (this->parallel_well_info_->communication().size() == 1)
                {
                    // Do the same thing as before. The else branch
                    // produces different rounding errors and results
                    // slightly different iteration counts / well curves
                    B_->mmv(x, y);
                }
                else
                {
                    Y temp(y);
                    mv(x, temp); // includes parallel reduction
                    y -= temp;
                }
            }
        private:
            const Matrix* B_;
            const ParallelWellInfo* parallel_well_info_;
        };

        inline
        double computeHydrostaticCorrection(const double well_ref_depth, const double vfp_ref_depth,
                                            const double rho, const double gravity) {
            const double dh = vfp_ref_depth - well_ref_depth;
            const double dp = rho * gravity * dh;

            return dp;
        }



        /// \brief Sums entries of the diagonal Matrix for distributed wells
        template<typename Scalar, typename Comm>
        void sumDistributedWellEntries(Dune::DynamicMatrix<Scalar>& mat, Dune::DynamicVector<Scalar>& vec,
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



        template <int dim, class C2F, class FC>
        std::array<double, dim>
        getCubeDim(const C2F& c2f,
                   FC         begin_face_centroids,
                   int        cell)
        {
            std::array< std::vector<double>, dim > X;
            {
                const std::vector<double>::size_type
                    nf = std::distance(c2f[cell].begin(),
                                       c2f[cell].end  ());

                for (int d = 0; d < dim; ++d) {
                    X[d].reserve(nf);
                }
            }

            typedef typename C2F::row_type::const_iterator FI;

            for (FI f = c2f[cell].begin(), e = c2f[cell].end(); f != e; ++f) {
                using Opm::UgGridHelpers::increment;
                using Opm::UgGridHelpers::getCoordinate;

                const FC& fc = increment(begin_face_centroids, *f, dim);

                for (int d = 0; d < dim; ++d) {
                    X[d].push_back(getCoordinate(fc, d));
                }
            }

            std::array<double, dim> cube;
            for (int d = 0; d < dim; ++d) {
                typedef std::vector<double>::iterator VI;
                typedef std::pair<VI,VI>              PVI;

                const PVI m = std::minmax_element(X[d].begin(), X[d].end());

                cube[d] = *m.second - *m.first;
            }

            return cube;
        }

        // explicite transpose of  dense matrix due to compilation problems
        // used for caclulating quasiimpes well weights
        template <class DenseMatrix>
        DenseMatrix transposeDenseDynMatrix(const DenseMatrix& M)
        {
            DenseMatrix tmp{M.cols(), M.rows()};
            for (size_t i = 0; i < M.rows(); ++i) {
                for (size_t j = 0; j < M.cols(); ++j) {
                    tmp[j][i] = M[i][j];
                }
            }
            return tmp;
        }

    } // namespace wellhelpers
}

#endif
