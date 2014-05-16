/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#include "config.h"
#include <opm/core/linalg/LinearSolverPetsc.hpp>
#include <opm/core/linalg/call_petsc.h>
#include <petsc.h>
#include <opm/core/utility/ErrorMacros.hpp>
namespace Opm
{

    LinearSolverPetsc::LinearSolverPetsc()
    {
        OPM_THROW(std::runtime_error, "Pestc just can be called through paramers.\n");
    }


    LinearSolverPetsc::LinearSolverPetsc(const parameter::ParameterGroup& param)
        : ksp_type_("gmres")
        , pc_type_("sor")
        , view_ksp_(false)
        , rtol_(1e-5)
        , atol_(1e-50)
        , dtol_(1e5)
        , maxits_(1e5)
    {
        int argc = 0;
        char** argv = NULL;
        PetscInitialize(&argc, &argv, (char*)0, "Petsc interface for OPM!\n");
        ksp_type_ = (param.getDefault("ksp_type", std::string(ksp_type_)));
        pc_type_ = (param.getDefault("pc_type", std::string(pc_type_)));
        view_ksp_ = (param.getDefault("ksp_view", int(view_ksp_)));
        rtol_ = param.getDefault("ksp_rtol", rtol_);
        atol_ = param.getDefault("ksp_atol", atol_);
        dtol_ = param.getDefault("ksp_dtol", dtol_);
        maxits_ = param.getDefault("ksp_max_it", maxits_);
    }


    LinearSolverPetsc::~LinearSolverPetsc()
    {
       PetscFinalize();
    }



    LinearSolverInterface::LinearSolverReport
    LinearSolverPetsc::solve(const int size,
                               const int nonzeros,
                               const int* ia,
                               const int* ja,
                               const double* sa,
                               const double* rhs,
                               double* solution,
                               const boost::any&) const
    {
        int ksp_type=4;
        int pc_type=0;
        const std::string ksp[]={"richardson", "chebyshev", "cg", "bicg", 
                                 "gmres", "fgmres", "dgmres", "gcr", "bcgs", 
                                 "cgs", "tfqmr", "tcqmr", "cr", "lsqr", "preonly"};
        const std::string pc[] = {"jacobi", "bjacobi", "sor", "eisenstat", 
                                  "icc", "ilu", "pcasm", "gamg", "ksp", 
                                  "composit", "lu", "cholesky", "none"};
        for (int i = 0; i < 15; ++i){
            if (ksp[i] == ksp_type_){
                ksp_type=i;
                break;        
            }
        }
        for (int i = 0; i < 13; ++i){
            if (pc[i] == pc_type_){
                pc_type=i;
                break;        
            }
        }
        call_Petsc(size, nonzeros, ia, ja, sa, rhs, solution, ksp_type, pc_type, rtol_, atol_, dtol_, maxits_, view_ksp_);
        LinearSolverReport rep = {};
        rep.converged = true;
        return rep;
    }

    void LinearSolverPetsc::setTolerance(const double /*tol*/)
    {
    }

    double LinearSolverPetsc::getTolerance() const
    {
        return -1.;
    }


} // namespace Opm

