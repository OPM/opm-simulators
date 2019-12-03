/*
  Copyright 2019 Big Data Accelerate

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
#include <memory>

#include <opm/bda/BdaBridge.hpp>
#include <opm/bda/BdaResult.hpp>

#define PRINT_TIMERS_BRIDGE_BRIDGE 0

typedef Dune::InverseOperatorResult InverseOperatorResult;

namespace Opm
{

BdaBridge::BdaBridge(bool use_gpu_, int maxit, double tolerance) : use_gpu(use_gpu_){
#if HAVE_CUDA
    if(use_gpu){
    	backend = new cusparseSolverBackend(maxit, tolerance);
    }
#endif
}

BdaBridge::~BdaBridge(){
#if HAVE_CUDA
	if(use_gpu){
		delete backend;
	}
#endif
}

#if HAVE_CUDA
template <class BridgeMatrix>
int checkZeroDiagonal(BridgeMatrix& mat) {
	static std::vector<int> diag_indices;   // contains offsets of the diagonal nnzs
	int numZeros = 0;
	const int dim = 3;
	const double zero_replace = 1e-15;
	double *nnzs = &(mat[0][0][0][0]);
	if(diag_indices.size() == 0){
		int N = mat.N()*dim;
		diag_indices.reserve(N);
		for(typename BridgeMatrix::const_iterator r = mat.begin(); r != mat.end(); ++r){
			for(auto c = r->begin(); c != r->end(); ++c){
				if(r.index() == c.index()){
					for(int rr = 0; rr < dim; ++rr){
						// pointer arithmetic
						int offset = (int)((long unsigned)&(mat[r.index()][c.index()][rr][rr]) - (long unsigned)nnzs); // in bytes
						offset /= sizeof(double);  // convert offset to doubles
						diag_indices.emplace_back(offset);
						double val = nnzs[offset];
						if(val == 0.0){ // could be replaced by '< 1e-30' or similar
							nnzs[offset] = zero_replace;
							++numZeros;
						}
					}
					break;
				}
			}
		}
	}else{
		for(int offset : diag_indices){
			if(nnzs[offset] == 0.0){ // could be replaced by '< 1e-30' or similar
				nnzs[offset] = zero_replace;
				++numZeros;
			}
		}
	}
	return numZeros;
}


// convert matrix to blocked csr (bsr) arrays
// if only_vals, do not convert rowPointers and colIndices
// sparsity pattern should stay the same due to matrix-add-well-contributions
template <class BridgeMatrix>
void convertMatrixBsr(BridgeMatrix& mat, std::vector<double> &h_vals, std::vector<int> &h_rows, std::vector<int> &h_cols, int dim, bool only_vals) {
	int sum_nnzs = 0;
	int nnz = mat.nonzeroes()*dim*dim;

	// copy nonzeros
	memcpy(h_vals.data(), &(mat[0][0][0][0]), sizeof(double)*nnz);

	// convert colIndices and rowPointers
	if(only_vals == false){
		h_rows.emplace_back(0);
			for(typename BridgeMatrix::const_iterator r = mat.begin(); r != mat.end(); ++r){
				int size_row = 0;
				for(auto c = r->begin(); c != r->end(); ++c){
					h_cols.emplace_back(c.index());
					size_row++;
				}
				sum_nnzs += size_row;
				h_rows.emplace_back(sum_nnzs);
			}
		// set last rowpointer
		h_rows[mat.N()] = mat.nonzeroes();
	}
} // end convertMatrixBsr()

// converts the BlockVector b to a flat array
template <class BridgeVector>
void convertBlockVectorToArray(BridgeVector& b, std::vector<double> &h_b) {
	memcpy(h_b.data(), &(b[0]), sizeof(double) * b.N() * b[0].dim());
}
#endif

template <class BridgeMatrix, class BridgeVector>
void BdaBridge::solve_system(BridgeMatrix *mat, BridgeVector &b, InverseOperatorResult &res)
{

#if HAVE_CUDA
	if(use_gpu){
		BdaResult result;
		result.converged = false;
		static std::vector<double> h_vals;
		static std::vector<double> h_b;
		static std::vector<int> h_rows;
		static std::vector<int> h_cols;
		int dim = (*mat)[0][0].N();
		int N = mat->N()*dim;
		int nnz = mat->nonzeroes()*dim*dim;

		if(dim != 3){
			std::cerr << "Error can only use cusparseSolver with blocksize = 3 at this time" << std::endl;
			exit(1);
		}

		if(h_vals.capacity() == 0){
			h_vals.reserve(nnz);
			h_vals.resize(nnz);
			h_b.reserve(N);
			h_b.resize(N);
			h_rows.reserve(N+1);
			h_cols.reserve(nnz);
		}

#if PRINT_TIMERS_BRIDGE
		Dune::Timer t_zeros;
		int numZeros = checkZeroDiagonal(*mat);
		printf("Checking zeros took %f s, found %d zeros\n", t_zeros.stop(), numZeros);
#else
		checkZeroDiagonal(*mat);
#endif

		bool initialized = backend->isInitialized();

#if PRINT_TIMERS_BRIDGE
		Dune::Timer t;
#endif

		convertMatrixBsr(*mat, h_vals, h_rows, h_cols, dim, initialized);
		convertBlockVectorToArray(b, h_b);

#if PRINT_TIMERS_BRIDGE
		printf("Conversion to flat arrays: %.4f s\n", t.stop());
#endif

		/////////////////////////
		// actually solve

		if(initialized == false){
			backend->initialize(N, nnz, dim);
			backend->copy_system_to_gpu(h_vals.data(), h_rows.data(), h_cols.data(), h_b.data());
			backend->analyse_matrix();
		}else{
			backend->update_system_on_gpu(h_vals.data(), h_b.data());
		}
		backend->reset_prec_on_gpu();
		if(backend->create_preconditioner()){
			backend->solve_system(result);
		}
		res.iterations = result.iterations;
		res.reduction = result.reduction;
		res.converged = result.converged;
		res.conv_rate = result.conv_rate;
		res.elapsed = result.elapsed;
	}else{
		res.converged = false;
	}
#endif // HAVE_CUDA
}


template <class BridgeVector>
void BdaBridge::get_result(BridgeVector &x){
#if HAVE_CUDA
	if(use_gpu){
		double *h_x = backend->post_process();
        // convert flat array to blockvector
        memcpy(&(x[0]), h_x, sizeof(double) * x.N() * x[0].dim());
	}
#endif
}

#if HAVE_CUDA
template void BdaBridge::solve_system< \
Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 2, 2>, std::allocator<Ewoms::MatrixBlock<double, 2, 2> > > , \
Dune::BlockVector<Dune::FieldVector<double, 2>, std::allocator<Dune::FieldVector<double, 2> > > > \
(Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 2, 2>, std::allocator<Ewoms::MatrixBlock<double, 2, 2> > > *mat, \
	Dune::BlockVector<Dune::FieldVector<double, 2>, std::allocator<Dune::FieldVector<double, 2> > > &b, \
	InverseOperatorResult &res);

template void BdaBridge::solve_system< \
Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 3, 3>, std::allocator<Ewoms::MatrixBlock<double, 3, 3> > > , \
Dune::BlockVector<Dune::FieldVector<double, 3>, std::allocator<Dune::FieldVector<double, 3> > > > \
(Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 3, 3>, std::allocator<Ewoms::MatrixBlock<double, 3, 3> > > *mat, \
	Dune::BlockVector<Dune::FieldVector<double, 3>, std::allocator<Dune::FieldVector<double, 3> > > &b, \
	InverseOperatorResult &res);

template void BdaBridge::solve_system< \
Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 4, 4>, std::allocator<Ewoms::MatrixBlock<double, 4, 4> > > , \
Dune::BlockVector<Dune::FieldVector<double, 4>, std::allocator<Dune::FieldVector<double, 4> > > > \
(Dune::BCRSMatrix<Ewoms::MatrixBlock<double, 4, 4>, std::allocator<Ewoms::MatrixBlock<double, 4, 4> > > *mat, \
	Dune::BlockVector<Dune::FieldVector<double, 4>, std::allocator<Dune::FieldVector<double, 4> > > &b, \
	InverseOperatorResult &res);


template void BdaBridge::get_result< \
Dune::BlockVector<Dune::FieldVector<double, 2>, std::allocator<Dune::FieldVector<double, 2> > > > \
(Dune::BlockVector<Dune::FieldVector<double, 2>, std::allocator<Dune::FieldVector<double, 2> > > &x);

template void BdaBridge::get_result< \
Dune::BlockVector<Dune::FieldVector<double, 3>, std::allocator<Dune::FieldVector<double, 3> > > > \
(Dune::BlockVector<Dune::FieldVector<double, 3>, std::allocator<Dune::FieldVector<double, 3> > > &x);

template void BdaBridge::get_result< \
Dune::BlockVector<Dune::FieldVector<double, 4>, std::allocator<Dune::FieldVector<double, 4> > > > \
(Dune::BlockVector<Dune::FieldVector<double, 4>, std::allocator<Dune::FieldVector<double, 4> > > &x);
#endif


}


