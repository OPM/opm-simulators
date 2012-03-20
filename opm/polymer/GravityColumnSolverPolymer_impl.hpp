/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/polymer/GravityColumnSolverPolymer.hpp>
#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{

    template <class Model>
    GravityColumnSolverPolymer<Model>::GravityColumnSolverPolymer(Model& model,
						    const UnstructuredGrid& grid,
						    const double tol,
						    const int maxit)
	: model_(model), grid_(grid), tol_(tol), maxit_(maxit)
    {
    }

    namespace {
	struct ZeroVec
	{
	    double operator[](int) const { return 0.0; }
	};
	struct StateWithZeroFlux
	{
	    StateWithZeroFlux(std::vector<double>& s, std::vector<double>& c) : sat(s), cpoly(c) {}
	    const ZeroVec& faceflux() const { return zv; }
	    const std::vector<double>& saturation() const { return sat; }
	    std::vector<double>& saturation() { return sat; }
	    const std::vector<double>& concentration() const { return cpoly; }
	    std::vector<double>& concentration() { return cpoly; }
	    ZeroVec zv;
	    std::vector<double>& sat;
	    std::vector<double>& cpoly;
	};
	struct Vecs
	{
	    Vecs(int sz) : sol(sz, 0.0) {}
	    const std::vector<double>& solution() const { return sol; }
	    std::vector<double>& writableSolution() { return sol; }
	    std::vector<double> sol;
	};
	struct JacSys
	{
	    JacSys(int sz) : v(sz) {}
	    const Vecs& vector() const { return v; }
	    Vecs& vector() { return v; }
	    Vecs v;
	    typedef std::vector<double> vector_type;
	};
        
        struct BandMatrixCoeff 
        {
            BandMatrixCoeff(int N, int ku, int kl) : N_(N), ku_(ku), kl_(kl), nrow_(2*kl + ku + 1) {
            }
            
            
            // compute the position where to store the coefficient of a matrix A_{i,j} (i,j=0,...,N-1)
            // in a array which is sent to the band matrix solver of LAPACK.
            int operator ()(int i, int j) const {
                return kl_ + ku_ + i - j + j*nrow_;
            }

            const int N_;
            const int ku_;
            const int kl_;
            const int nrow_;
        };

    } // anon namespace



    /// \param[in] columns         for each column (with logical cartesian indices as key),
    ///                            contains the cells on which to solve the segregation
    ///                            problem. For each column, its cells must be in a single
    ///                            vertical column, and ordered
    ///                            (direction doesn't matter).
    template <class Model>
    void GravityColumnSolverPolymer<Model>::solve(const std::map<int, std::vector<int> >& columns,
						  const double dt,
						  std::vector<double>& s,
						  std::vector<double>& c
						  )
    {
	// Initialize model. These things are done for the whole grid!
	StateWithZeroFlux state(s, c); // This holds s and c by reference.
	JacSys sys(2*grid_.number_of_cells);
	std::vector<double> increment(2*grid_.number_of_cells, 0.0);
	model_.initStep(state, grid_, sys);

	int iter = 0;
	double max_delta = 1e100;
	while (iter < maxit_) {
	    model_.initIteration(state, grid_, sys);
	    std::map<int, std::vector<int> >::const_iterator it;
	    for (it = columns.begin(); it != columns.end(); ++it) {
		solveSingleColumn(it->second, dt, s, c, increment);
	    }
	    for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
		sys.vector().writableSolution()[2*cell + 0] += increment[2*cell + 0];
		sys.vector().writableSolution()[2*cell + 1] += increment[2*cell + 1];
	    }
	    const double maxelem = *std::max_element(increment.begin(), increment.end());
	    const double minelem = *std::min_element(increment.begin(), increment.end());
	    max_delta = std::max(maxelem, -minelem);
	    std::cout << "Iteration " << iter << "   max_delta = " << max_delta << std::endl;
	    if (max_delta < tol_) {
		break;
	    }
	    ++iter;
	}
	if (max_delta >= tol_) {
	    THROW("Failed to converge!");
	}
	// Finalize.
	// model_.finishIteration(); // Doesn't do anything in th 2p model.
	// finishStep() writes to state, which holds s by reference.
	// This will update the entire grid's state...
	model_.finishStep(grid_, sys.vector().solution(), state);
    }




    /// \param[in] column_cells    the cells on which to solve the segregation
    ///                            problem. Must be in a single vertical column,
    ///                            and ordered (direction doesn't matter).
    template <class Model>
    void GravityColumnSolverPolymer<Model>::solveSingleColumn(const std::vector<int>& column_cells,
                                                              const double dt,
                                                              std::vector<double>& s,
                                                              std::vector<double>& c,
                                                              std::vector<double>& sol_vec)
    {
	// This is written only to work with SinglePointUpwindTwoPhase,
	// not with arbitrary problem models.
	const int col_size = column_cells.size();
	StateWithZeroFlux state(s, c); // This holds s by reference.

	// Assemble.
        const int kl = 3;
        const int ku = 3;
        const int nrow = 2*kl + ku + 1;
        const int N = 2*col_size; // N unknowns: s and c for each cell.
        const BandMatrixCoeff bmc(N, ku, kl);

	std::vector<double> hm(nrow*N, 0.0); // band matrix with 3 upper and 3 lower diagonals.
	std::vector<double> rhs(N, 0.0);
	for (int ci = 0; ci < col_size; ++ci) {
	    std::vector<double> F(2, 0.);   
	    std::vector<double> dFd1(4, 0.);
	    std::vector<double> dFd2(4, 0.);
	    std::vector<double> dF(4, 0.);
	    const int cell = column_cells[ci];
	    const int prev_cell = (ci == 0) ? -999 : column_cells[ci - 1];
	    const int next_cell = (ci == col_size - 1) ? -999 : column_cells[ci + 1];
	    // model_.initResidual(cell, F);
	    for (int j = grid_.cell_facepos[cell]; j < grid_.cell_facepos[cell+1]; ++j) {
		const int face = grid_.cell_faces[j];
		const int c1 = grid_.face_cells[2*face + 0];
                const int c2 = grid_.face_cells[2*face + 1];
		if (c1 == prev_cell || c2 == prev_cell || c1 == next_cell || c2 == next_cell) {
                    F.assign(2, 0.);   
                    dFd1.assign(4, 0.);
                    dFd2.assign(4, 0.);
		    model_.fluxConnection(state, grid_, dt, cell, face, &F[0], &dFd1[0], &dFd2[0]);
		    if (c1 == prev_cell || c2 == prev_cell) {
                        hm[bmc(2*ci + 0, 2*(ci - 1) + 0)] += dFd2[0];
                        hm[bmc(2*ci + 0, 2*(ci - 1) + 1)] += dFd2[1];
                        hm[bmc(2*ci + 1, 2*(ci - 1) + 0)] += dFd2[2];
                        hm[bmc(2*ci + 1, 2*(ci - 1) + 1)] += dFd2[3];
		    } else {
			ASSERT(c1 == next_cell || c2 == next_cell);
                        hm[bmc(2*ci + 0, 2*(ci + 1) + 0)] += dFd2[0];
                        hm[bmc(2*ci + 0, 2*(ci + 1) + 1)] += dFd2[1];
                        hm[bmc(2*ci + 1, 2*(ci + 1) + 0)] += dFd2[2];
                        hm[bmc(2*ci + 1, 2*(ci + 1) + 1)] += dFd2[3];
		    }
                    hm[bmc(2*ci + 0, 2*ci + 0)] += dFd1[0];
                    hm[bmc(2*ci + 0, 2*ci + 1)] += dFd1[1];
                    hm[bmc(2*ci + 1, 2*ci + 0)] += dFd1[2];
                    hm[bmc(2*ci + 1, 2*ci + 1)] += dFd1[3];
                    
		    rhs[2*ci + 0] += F[0];
		    rhs[2*ci + 1] += F[1];
		}
	    }
	    F.assign(2, 0.);
            dF.assign(4, 0.);
	    model_.accumulation(grid_, cell, &F[0], &dF[0]);
            hm[bmc(2*ci + 0, 2*ci + 0)] += dF[0];
            hm[bmc(2*ci + 0, 2*ci + 1)] += dF[1];
            hm[bmc(2*ci + 1, 2*ci + 0)] += dF[2];
            hm[bmc(2*ci + 1, 2*ci + 1)] += dF[3];

            rhs[2*ci + 0] += F[0];
            rhs[2*ci + 1] += F[1];

	}
	// model_.sourceTerms(); // Not needed 
	// Solve.
	const int num_rhs = 1;
	int info = 0;
        int ipiv;
	// Solution will be written to rhs.
        dgbsv_(&N, &kl, &ku, &num_rhs, &hm[0], &nrow, &ipiv, &rhs[0], &N, &info);
	if (info != 0) {
	    THROW("Lapack reported error in dgtsv: " << info);
	}
	for (int ci = 0; ci < col_size; ++ci) {
	    sol_vec[2*column_cells[ci] + 0] = -rhs[2*ci + 0];
	    sol_vec[2*column_cells[ci] + 1] = -rhs[2*ci + 1];
	}
    }

} // namespace Opm
