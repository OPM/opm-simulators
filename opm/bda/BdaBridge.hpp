
#ifndef BDABRIDGE_HEADER_INCLUDED
#define BDABRIDGE_HEADER_INCLUDED

#include "dune/istl/solver.hh" // for struct InverseOperatorResult

#include "dune/istl/bcrsmatrix.hh"
#include <ewoms/linear/matrixblock.hh>

#ifdef __NVCC__
#include <opm/bda/cusparseSolverBackend.hpp>
#endif

namespace Opm
{

typedef Dune::InverseOperatorResult InverseOperatorResult;


class BdaBridge
{
private:
#ifdef __NVCC__
	cusparseSolverBackend *backend;
#endif
	bool use_gpu;

public:
	BdaBridge(bool use_gpu, int maxit, double tolerance);

	~BdaBridge();

	template <class BridgeMatrix, class BridgeVector>
	void solve_system(BridgeMatrix *mat, BridgeVector &b, InverseOperatorResult &result);

	template <class BridgeVector>
	void get_result(BridgeVector &x);

}; // end class BdaBridge

}

#endif
