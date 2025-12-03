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
#ifndef ELASTICITY_PRIMARY_VARIABLES_HPP
#define ELASTICITY_PRIMARY_VARIABLES_HPP

#include <dune/common/fvector.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/utils/basicproperties.hh>

#include <type_traits>


namespace Opm {

template <class TypeTag>
class ElasticityPrimaryVariables
    : public Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                               getPropValue<TypeTag, Properties::NumEqTPSA>()>
{
    using Evaluation = GetPropType<TypeTag, Properties::EvaluationTPSA>;
    using Indices = GetPropType<TypeTag, Properties::IndicesTPSA>;
    using Scalar = GetPropType<TypeTag, Opm::Properties::Scalar>;

    using Toolbox = Opm::MathToolbox<Evaluation>;

    enum { disp0Idx = Indices::disp0Idx };
    enum { rot0Idx = Indices::rot0Idx };
    enum { solidPres0Idx = Indices::solidPres0Idx };

    enum { numEq = getPropValue<TypeTag, Properties::NumEqTPSA>() };

    using ParentType = Dune::FieldVector<Scalar, numEq>;

public:
    /*!
    * \brief Constructor
    */
    ElasticityPrimaryVariables() : ParentType()
    {
        Valgrind::SetUndefined(*this);
    }

    /*!
    * \brief Default copy constructor
    */
    ElasticityPrimaryVariables(const ElasticityPrimaryVariables& value) = default;

    /*!
    * \brief Default assignment constructor
    */
    ElasticityPrimaryVariables& operator=(const ElasticityPrimaryVariables& value) = default;

    using ParentType::operator=; //!< Import base class assignment operators.

    /*!
    * \brief Return primary variable in Evaluation type
    *
    * \param varIdx Primary variable index
    * \param timeIdx Time index
    * \param linearizationType Type of linearization
    *
    * Automatic differentiation:    returns value + derivative
    * Finite differences:           returns value only
    */
    Evaluation makeEvaluation(unsigned varIdx,
                              unsigned timeIdx,
                              Opm::LinearizationType linearizationType = LinearizationType()) const
    {
        // Finite difference
        if constexpr (std::is_same_v<Evaluation, Scalar>) {
            return (*this)[varIdx];
        }
        // Automatic differentiation
        else {
            if (timeIdx == linearizationType.time) {
                return Toolbox::createVariable((*this)[varIdx], varIdx);
            }
            else {
                return Toolbox::createConstant((*this)[varIdx]);
            }
        }
    }

    /*!
    * \brief Assign primary variables from a material state container
    *
    * \param materialState Material state container
    */
    template <class MaterialState>
    void assignNaive(const MaterialState& materialState)
    {
        // Reset primary variables vector
        (*this) = 0.0;

        // Assign displacement and rotation vectors
        for (unsigned dirIdx = 0; dirIdx < 3; ++dirIdx) {
            (*this)[disp0Idx + dirIdx] = Opm::getValue(materialState.displacement(dirIdx));
            (*this)[rot0Idx + dirIdx] = Opm::getValue(materialState.rotation(dirIdx));
        }

        // Assign solid pressure
        (*this)[solidPres0Idx] = Opm::getValue(materialState.solidPressure());
    }

    /*!
    * \brief Instruct Valgrind to check the definedness of all attributes of this class
    */
    void checkDefined() const
    {
        Valgrind::CheckDefined(*static_cast<const ParentType*>(this));
    }
};  // class ElasticityPrimaryVariables

}  // namespace Opm

#endif
