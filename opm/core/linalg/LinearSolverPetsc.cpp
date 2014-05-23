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
#include <unordered_map>
#include <petsc.h>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{
    
namespace{

    class KSPTypeMap {
    public:
        explicit
        KSPTypeMap(const std::string& default_type = "gmres")
            : default_type_(default_type)
        {
            type_map_.insert(std::make_pair("richardson", KSPRICHARDSON));
            type_map_.insert(std::make_pair("chebyshev", KSPCHEBYSHEV));
            type_map_.insert(std::make_pair("cg", KSPCG));
            type_map_.insert(std::make_pair("bicgs", KSPBICG));
            type_map_.insert(std::make_pair("gmres", KSPGMRES));
            type_map_.insert(std::make_pair("fgmres", KSPFGMRES));
            type_map_.insert(std::make_pair("dgmres", KSPDGMRES));
            type_map_.insert(std::make_pair("gcr", KSPGCR));
            type_map_.insert(std::make_pair("bcgs", KSPBCGS));
            type_map_.insert(std::make_pair("cgs", KSPCGS));
            type_map_.insert(std::make_pair("tfqmr", KSPTFQMR));
            type_map_.insert(std::make_pair("tcqmr", KSPTCQMR));
            type_map_.insert(std::make_pair("cr", KSPCR));
            type_map_.insert(std::make_pair("preonly", KSPPREONLY));
        }
          
        KSPType 
        find(const std::string& type) const
        {
            Map::const_iterator it = type_map_.find(type);

            if (it == type_map_.end()) {
                it = type_map_.find(default_type_);
            }
        
            if (it == type_map_.end()) {
                OPM_THROW(std::runtime_error, "Unknown KSPType: '" << type << "'");
            }
            return it->second;
        }
    private:
        typedef std::unordered_map<std::string, KSPType> Map;

        std::string default_type_;
        Map type_map_;
    };


    class PCTypeMap {
    public:
        explicit
        PCTypeMap(const std::string& default_type = "jacobi")
            : default_type_(default_type)
        {
            type_map_.insert(std::make_pair("jacobi", PCJACOBI));
            type_map_.insert(std::make_pair("bjacobi", PCBJACOBI));
            type_map_.insert(std::make_pair("sor", PCSOR));
            type_map_.insert(std::make_pair("eisenstat", PCEISENSTAT));
            type_map_.insert(std::make_pair("icc", PCICC));
            type_map_.insert(std::make_pair("ilu", PCILU));
            type_map_.insert(std::make_pair("asm", PCASM));
            type_map_.insert(std::make_pair("gamg", PCGAMG));
            type_map_.insert(std::make_pair("ksp", PCKSP));
            type_map_.insert(std::make_pair("composite", PCCOMPOSITE));
            type_map_.insert(std::make_pair("lu", PCLU));
            type_map_.insert(std::make_pair("cholesky", PCCHOLESKY));
            type_map_.insert(std::make_pair("none", PCNONE));
        }
            
        PCType 
        find(const std::string& type) const
        {
            Map::const_iterator it = type_map_.find(type);

            if (it == type_map_.end()) {
                it = type_map_.find(default_type_);
            }
       
            if (it == type_map_.end()) {
                OPM_THROW(std::runtime_error, "Unknown PCType: '" << type << "'");
            }
            return it->second;
        }
    private:
        typedef std::unordered_map<std::string, PCType> Map;

        std::string default_type_;
        Map type_map_;
    };

} // anonymous namespace.

    LinearSolverPetsc::LinearSolverPetsc()
    {
        OPM_THROW(std::runtime_error, "Pestc just can be called through paramers.\n");
    }


    LinearSolverPetsc::LinearSolverPetsc(const parameter::ParameterGroup& param)
        : ksp_type_("gmres")
        , pc_type_("sor")
        , ksp_view_(false)
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
        ksp_view_ = (param.getDefault("ksp_view", int(ksp_view_)));
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
        KSPTypeMap ksp(ksp_type_);
        KSPType ksp_type = ksp.find(ksp_type_);
        PCTypeMap pc(pc_type_);
        PCType pc_type = pc.find(pc_type_);
        
        call_Petsc(size, nonzeros, ia, ja, sa, rhs, solution, ksp_type, pc_type, rtol_, atol_, dtol_, maxits_, ksp_view_);

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

