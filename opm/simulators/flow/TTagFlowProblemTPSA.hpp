// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef TTAG_FLOW_PROBLEM_TPSA_HPP
#define TTAG_FLOW_PROBLEM_TPSA_HPP

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/material/densead/Evaluation.hpp>

#include <opm/models/discretization/common/tpsalinearizer.hpp>
#include <opm/models/tpsa/elasticityindices.hpp>
#include <opm/models/tpsa/elasticitylocalresidualtpsa.hpp>
#include <opm/models/tpsa/elasticityprimaryvariables.hpp>
#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/tpsa/tpsamodel.hpp>
#include <opm/models/tpsa/tpsanewtonconvergencewriter.hpp>
#include <opm/models/tpsa/tpsanewtonmethod.hpp>

#include <opm/simulators/flow/BlackoilModelTPSA.hpp>
#include <opm/simulators/flow/FlowProblemTPSA.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/ISTLSolverTPSA.hpp>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>


namespace Opm::Properties {

namespace TTag {

struct FlowProblemTpsa {};

}  // Opm::Properties::TTag

// TPSA indices for primary variables and equations
template<class TypeTag>
struct IndicesTPSA<TypeTag, TTag::FlowProblemTpsa>
{
    using type = ElasticityIndices</*PVOffset=*/0>;
};

// Number of TPSA equations
template<class TypeTag>
struct NumEqTPSA<TypeTag, TTag::FlowProblemTpsa>
{ static constexpr int value = GetPropType<TypeTag, Properties::IndicesTPSA>::numEq; };

// TPSA linearizer
template<class TypeTag>
struct LinearizerTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = TpsaLinearizer<TypeTag>; };

// Set the function evaluation w.r.t. the TPSA primary variables
template<class TypeTag>
struct EvaluationTPSA<TypeTag, TTag::FlowProblemTpsa>
{
private:
    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEqTPSA>();

    using Scalar = GetPropType<TypeTag, Scalar>;

public:
    using type = DenseAd::Evaluation<Scalar, numEq>;
};

// TPSA equation vector
template<class TypeTag>
struct EqVectorTPSA<TypeTag, TTag::FlowProblemTpsa>
{
    using type = Dune::FieldVector<GetPropType<TypeTag, Scalar>,
                                   getPropValue<TypeTag, Properties::NumEqTPSA>()>;
};

// Global TPSA equation vector
template<class TypeTag>
struct GlobalEqVectorTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::EqVectorTPSA>>; };

// TPSA Newton method
template<class TypeTag>
struct NewtonMethodTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = TpsaNewtonMethod<TypeTag>; };

template<class TypeTag>
struct NewtonConvergenceWriterTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = TpsaNewtonConvergenceWriter<TypeTag>; };

// TPSA primary variables
template<class TypeTag>
struct PrimaryVariablesTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = ElasticityPrimaryVariables<TypeTag>; };

// TPSA solution vector
template<class TypeTag>
struct SolutionVectorTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariablesTPSA>>; };

// TPSA number of historic solutions to save
template<class TypeTag>
struct SolutionHistorySizeTPSA<TypeTag, TTag::FlowProblemTpsa>
{ static constexpr int value = 2; };

// TPSA model
template<class TypeTag>
struct ModelTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = TpsaModel<TypeTag>; };

// TPSA local residual
template<class TypeTag>
struct LocalResidualTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = ElasticityLocalResidual<TypeTag>; };

// TPSA sparse matrix adapter for Jacobian
template<class TypeTag>
struct SparseMatrixAdapterTPSA<TypeTag, TTag::FlowProblemTpsa>
{
private:
    using Scalar = GetPropType<TypeTag, Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;

};

// Disable constraints in Newton method
template<class TypeTag>
struct EnableConstraintsTPSA<TypeTag, TTag::FlowProblemTpsa>
{ static constexpr bool value = false; };

// Set linear solver backend
template<class TypeTag>
struct LinearSolverBackendTPSA<TypeTag, TTag::FlowProblemTpsa>
{ using type = ISTLSolverTPSA<TypeTag>; };

}  // namespace Opm::Properties

#endif
