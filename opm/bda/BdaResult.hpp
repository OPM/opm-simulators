
#ifndef BDARESULT_HEADER_INCLUDED
#define BDARESULT_HEADER_INCLUDED

namespace Opm
{

// based on InverseOperatorResult struct from dune/istl/solver.hh
class BdaResult
{

public:
	int iterations = 0;			// number of iterations
	double reduction = 0.0;		// reduction of norm, norm_start / norm_final
	bool converged = false;		// true iff the linear solver reached the desired norm within maxit iterations
	double conv_rate = 0.0;		// average reduction of norm per iteration, usually calculated with 'static_cast<double>(pow(res.reduction,1.0/it));'
	double elapsed = 0.0;		// time in seconds to run the linear solver

	// Dune 2.6 has a member 'double condition_estimate = -1' in InverseOperatorResult

}; // end class BdaResult

}

#endif
