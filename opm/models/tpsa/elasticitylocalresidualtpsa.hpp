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
#ifndef ELASTICITY_LOCAL_RESIDUAL_TPSA_HPP
#define ELASTICITY_LOCAL_RESIDUAL_TPSA_HPP

#include <dune/common/fvector.hh>

#include <opm/input/eclipse/Schedule/BCProp.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/materialstates/MaterialStateTPSA.hpp>

#include <opm/models/tpsa/tpsabaseproperties.hpp>

#include <opm/simulators/flow/FacePropertiesTPSA.hpp>

#include <stdexcept>


namespace Opm {

template <class TypeTag>
class ElasticityLocalResidual
{
    using Indices = GetPropType<TypeTag, Properties::IndicesTPSA>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using MaterialState = MaterialStateTPSA<Evaluation>;
    using Toolbox = MathToolbox<Evaluation>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiRotEqIdx = Indices::contiRotEqIdx };
    enum { contiSolidPresEqIdx = Indices::contiSolidPresEqIdx };
    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };

public:
    /*!
    * \brief Calculate volume terms in TPSA formulation
    *
    * \param volTerm Volume term vector
    * \param materialState Material state container
    * \param problem Flow problem
    * \param globalIndex Cell index
    *
    * \note Material state, problem input and global index here might/should be merged in an "IntensiveQuantity"
    *       container as in BlackOilLocalResidualTPFA
    */
    template <class LhsEval>
    static void computeVolumeTerm(Dune::FieldVector<LhsEval, numEq>& volTerm,
                                  const MaterialState& materialState,
                                  const Problem& problem,
                                  const unsigned globalIndex)
    {
        // Reset volume terms
        volTerm = 0.0;

        // Rotation equations (one per direction)
        const Scalar sModulus = problem.shearModulus(globalIndex);
        for (unsigned dirIdx = 0; dirIdx < 3; ++dirIdx) {
            volTerm[contiRotEqIdx + dirIdx] +=
                Toolbox::template decay<LhsEval>(materialState.rotation(dirIdx))
                / sModulus;
        }

        // Solid pressure equation
        const Scalar lame = problem.lame(globalIndex);
        volTerm[contiSolidPresEqIdx] +=
            Toolbox::template decay<LhsEval>(materialState.solidPressure())
            / lame;
    }

    /*!
    * \brief Calculate terms across cell faces in TPSA formulation
    *
    * \param faceTerm Face term vector
    * \param materialStateIn Material state container of inside cell
    * \param materialStateEx Material state container of outside cel
    * \param problem Flow problem
    * \param globalIndexIn Inside cell index
    * \param globalIndexEx Outside cell index
    *
    * \note Material state, problem input and global index here might/should be merged in "IntensiveQuantity" and
    *       "NeighborInfo" containers as in BlackOilLocalResidualTPFA
    */
    static void computeFaceTerm(Dune::FieldVector<Evaluation, numEq>& faceTerm,
                                const MaterialState& materialStateIn,
                                const MaterialState& materialStateEx,
                                Problem& problem,
                                const unsigned globalIndexIn,
                                const unsigned globalIndexEx)
    {
        // Reset face terms
        faceTerm = 0.0;

        // Extract some face properties
        const Scalar weightAvgIn = problem.weightAverage(globalIndexIn, globalIndexEx);
        const Scalar weightAvgEx = problem.weightAverage(globalIndexEx, globalIndexIn);
        const Scalar weightProd = problem.weightProduct(globalIndexIn, globalIndexEx);
        const Scalar normDist = problem.normalDistance(globalIndexIn, globalIndexEx);
        const auto& faceNormal = problem.cellFaceNormal(globalIndexIn, globalIndexEx);

        // Effective shear modulus
        const Scalar sModulusIn = problem.shearModulus(globalIndexIn);
        const Scalar sModulusEx = problem.shearModulus(globalIndexEx);
        const Scalar eff_sModulus = weightAvgIn * sModulusIn + weightAvgEx * sModulusEx;

        // Distance ratio
        const Scalar distRatio = 0.5 * weightProd / normDist;

        // Solid pressures
        const Evaluation& solidPIn = materialStateIn.solidPressure();
        const Scalar solidPEx = decay<Scalar>(materialStateEx.solidPressure());

        // ///
        // Solid pressure equation (direction-independent equation)
        // ///
        faceTerm[contiSolidPresEqIdx] +=
            distRatio * eff_sModulus * (solidPIn - solidPEx);

        // ///
        // Displacement, rotation and solid pressure (directional-dependent) equations
        // ///
        // Lambda function for computing modulo 3 of possibly negative integers to get indices in cross product.
        // E.g. if i = x-dir(=0), we want y-dir(=1) and z-dir(=2), hence -1 mod 3 must equal 2 and not -1
        auto modNeg = [](int i) { return ((i % 3) + 3) % 3; };

        // Loop over x-, y- and z-dir (corresponding to dirIdx = 0, 1, 2)
        for (int dirIdx = 0; dirIdx < 3; ++dirIdx) {
            // Direction indices in cross-product
            unsigned dirIdxNeg = modNeg(dirIdx - 1);
            unsigned dirIdxPos = modNeg(dirIdx + 1);

            // Displacement equation
            const Scalar faceNormalDir = faceNormal[dirIdx];
            const Scalar faceNormalNeg = faceNormal[dirIdxNeg];
            const Scalar faceNormalPos = faceNormal[dirIdxPos];

            const Evaluation& dispIn = materialStateIn.displacement(dirIdx);
            const Scalar dispEx = decay<Scalar>(materialStateEx.displacement(dirIdx));

            const Evaluation& rotInNeg = materialStateIn.rotation(dirIdxNeg);
            const Evaluation& rotInPos = materialStateIn.rotation(dirIdxPos);
            const Scalar rotExNeg =  decay<Scalar>(materialStateEx.rotation(dirIdxNeg));
            const Scalar rotExPos =  decay<Scalar>(materialStateEx.rotation(dirIdxPos));

            faceTerm[conti0EqIdx + dirIdx] +=
                2.0 * (eff_sModulus / normDist) * (dispIn - dispEx)
                - weightAvgIn * (faceNormalNeg * rotInPos - faceNormalPos * rotInNeg)
                - weightAvgEx * (faceNormalNeg * rotExPos - faceNormalPos * rotExNeg)
                - faceNormalDir * (weightAvgIn * solidPIn + weightAvgEx * solidPEx);

            // Rotation equation
            const Evaluation& dispInNeg = materialStateIn.displacement(dirIdxNeg);
            const Evaluation& dispInPos = materialStateIn.displacement(dirIdxPos);
            const Scalar dispExNeg = decay<Scalar>(materialStateEx.displacement(dirIdxNeg));
            const Scalar dispExPos = decay<Scalar>(materialStateEx.displacement(dirIdxPos));

            faceTerm[contiRotEqIdx + dirIdx] +=
                - weightAvgEx * (faceNormalNeg * dispInPos - faceNormalPos * dispInNeg)
                - weightAvgIn * (faceNormalNeg * dispExPos - faceNormalPos * dispExNeg);

            // Solid pressure equation
            faceTerm[contiSolidPresEqIdx] +=
                - faceNormalDir * (weightAvgEx * dispIn + weightAvgIn * dispEx);
        }
    }

    /*!
    * \brief Calculate boundary conditions in TPSA formulation given by BCCON/BCPROP
    *
    * \param bndryTerm Boundary term vector
    * \param materialState Material state container
    * \param bdyInfo Boundary condition info container
    * \param problem Flow problem
    * \param globalIndex Cell index
    */
    template <class BoundaryConditionData>
    static void computeBoundaryTerm(Dune::FieldVector<Evaluation, numEq>& bndryTerm,
                                    const MaterialState& materialState,
                                    const BoundaryConditionData& bdyInfo,
                                    Problem& problem,
                                    unsigned globalIndex)
    {
        // Switch between possible boundary conditions
        switch (bdyInfo.type) {
        // OBS: NONE is interpreted as FIXED with zero displacement
        case BCMECHType::NONE:
            computeBoundaryTermFixed(bndryTerm,
                                     materialState,
                                     bdyInfo,
                                     problem,
                                     globalIndex);
            break;
        case BCMECHType::FIXED:
            throw std::runtime_error("BCTYPE FIXED has not been implemented in TPSA");
        case BCMECHType::FREE:
            computeBoundaryTermFree(bndryTerm,
                                    materialState,
                                    bdyInfo,
                                    problem,
                                    globalIndex);
            break;
        default:
            throw std::logic_error("Unknown boundary condition type " +
                                    std::to_string(static_cast<int>(bdyInfo.type)) +
                                    " in computeBoundaryFlux()." );
        }
    }

    /*!
    * \brief Calculate fixed displacement boundary condition in TPSA formulation
    *
    * \param bndryTerm Boundary term vector
    * \param materialState Material state container
    * \param bdyInfo Boundary condition info container
    * \param problem Flow problem
    * \param globalIndex Cell index
    *
    * \note BCMECHType::NONE is a implemented as a specialization of BCMECHType::FIXED with displacement equal
    * to zero on the boundary face
    */
    template <class BoundaryConditionData>
    static void computeBoundaryTermFixed(Dune::FieldVector<Evaluation, numEq>& bndryTerm,
                                         const MaterialState& materialState,
                                         const BoundaryConditionData& bdyInfo,
                                         Problem& problem,
                                         unsigned globalIndex)
    {
        // !!!!
        // Only BCMECHType::NONE, where we have zero displacement on boundary face, have been implemented!
        // !!!!

        // Reset bondary term
        bndryTerm = 0.0;

        // Extract some face properties
        const unsigned bfIdx = bdyInfo.boundaryFaceIndex;
        const Scalar weightAvg = problem.weightAverageBoundary(globalIndex, bfIdx);
        const Scalar normDist = problem.normalDistanceBoundary(globalIndex, bfIdx);
        const auto& faceNormal = problem.cellFaceNormalBoundary(globalIndex, bfIdx);

        // Effective shear modulus (= cell shear modulus)
        const Scalar eff_sModulus = problem.shearModulus(globalIndex);

        // Solid pressure
        const Evaluation& solidP = materialState.solidPressure();

        // ///
        // Displacement equation
        // ///
        // Lambda function for computing modulo 3 of possibly negative integers to get indices in cross product.
        // E.g. if i = x-dir(=0), we want y-dir(=1) and z-dir(=2), hence -1 mod 3 must equal 2 and not -1
        auto modNeg = [](int i) { return ((i % 3) + 3) % 3; };

        // Loop over x-, y- and z-dir (corresponding to dirIdx = 0, 1, 2)
        for (int dirIdx = 0; dirIdx < 3; ++dirIdx) {
            // Direction indices in cross-product
            unsigned dirIdxNeg = modNeg(dirIdx - 1);
            unsigned dirIdxPos = modNeg(dirIdx + 1);

            // Displacement equation
            const Scalar faceNormalDir = faceNormal[dirIdx];
            const Scalar faceNormalNeg = faceNormal[dirIdxNeg];
            const Scalar faceNormalPos = faceNormal[dirIdxPos];

            const Evaluation& disp = materialState.displacement(dirIdx);

            const Evaluation& rotNeg = materialState.rotation(dirIdxNeg);
            const Evaluation& rotPos = materialState.rotation(dirIdxPos);

            bndryTerm[conti0EqIdx + dirIdx] +=
                2.0 * (eff_sModulus / normDist) * disp
                - weightAvg * (faceNormalNeg * rotPos - faceNormalPos * rotNeg)
                - faceNormalDir * weightAvg * solidP;
        }
    }

    /*!
    * \brief Calculate free (or zero traction) boundary condition in TPSA formulation
    *
    * \param bndryTerm Boundary term vector
    * \param materialState Material state container
    * \param bdyInfo Boundary condition info container
    * \param problem Flow problem
    * \param globalIndex Cell index
    *
    * \note Free, or zero traction, BC is equivalent of having a spring at infinity where we have assumed all (primary)
    * variables and parameters (e.g., shear modulus) are zero.
    */
    template <class BoundaryConditionData>
    static void computeBoundaryTermFree(Dune::FieldVector<Evaluation, numEq>& bndryTerm,
                                        const MaterialState& materialState,
                                        const BoundaryConditionData& bdyInfo,
                                        Problem& problem,
                                        unsigned globalIndex)
    {
        // Reset bondary term
        bndryTerm = 0.0;

        // Face properties
        const unsigned bfIdx = bdyInfo.boundaryFaceIndex;
        const Scalar weightAvg = 1.0;
        const Scalar normDist = problem.normalDistanceBoundary(globalIndex, bfIdx);
        const auto& faceNormal = problem.cellFaceNormalBoundary(globalIndex, bfIdx);

        const Scalar sModulus = problem.shearModulus(globalIndex);

        // ///
        // Solid pressure equation (direction-independent equation)
        // ///
        const Evaluation& solidP = materialState.solidPressure();
        bndryTerm[contiSolidPresEqIdx] +=
            0.5 * (normDist / sModulus) * solidP;

        // ///
        // Rotation and solid pressure (directional-dependent) equations
        // ///
        // Lambda function for computing modulo 3 of possibly negative integers to get indices in cross product.
        // E.g. if i = x-dir(=0), we want y-dir(=1) and z-dir(=2), hence -1 mod 3 must equal 2 and not -1
        auto modNeg = [](int i) { return ((i % 3) + 3) % 3; };

        // Loop over x-, y- and z-dir (corresponding to dirIdx = 0, 1, 2)
        for (int dirIdx = 0; dirIdx < 3; ++dirIdx) {
            // Direction indices in cross-product
            unsigned dirIdxNeg = modNeg(dirIdx - 1);
            unsigned dirIdxPos = modNeg(dirIdx + 1);

            // Rotation equation
            const Scalar faceNormalNeg = faceNormal[dirIdxNeg];
            const Scalar faceNormalPos = faceNormal[dirIdxPos];

            const Evaluation& dispNeg = materialState.displacement(dirIdxNeg);
            const Evaluation& dispPos = materialState.displacement(dirIdxPos);

            bndryTerm[contiRotEqIdx + dirIdx] +=
                - weightAvg * (faceNormalNeg * dispPos - faceNormalPos * dispNeg);

            // Solid pressure (directional-dependent) equation
            const Scalar faceNormalDir = faceNormal[dirIdx];
            const Evaluation& disp = materialState.displacement(dirIdx);

            bndryTerm[contiSolidPresEqIdx] +=
                - faceNormalDir * weightAvg * disp;
        }
    }

    /*!
    * \brief Calculate source term in TPSA formulation
    *
    * \param sourceTerm Source term vector
    * \param problem Flow problem
    * \param globalIndex Cell index
    * \param timeIdx Time index
    */
    static void computeSourceTerm(Dune::FieldVector<Evaluation, numEq>& sourceTerm,
                                  Problem& problem,
                                  unsigned globalSpaceIdex,
                                  unsigned timeIdx)
    {
        // Reset source terms
        sourceTerm = 0.0;

        // Get source term from problem
        // NOTE: separate source function than Flow source(...)!
        problem.tpsaSource(sourceTerm, globalSpaceIdex, timeIdx);
    }
};  // class ElasticityLocalResidual

}  // namespace Opm

#endif