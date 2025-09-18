/*
  Copyright 2024,2025 Equinor ASA.

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

#ifndef OPM_INSTANTIATION_INDICES_MACROS_HPP
#define OPM_INSTANTIATION_INDICES_MACROS_HPP

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

/**
 * \brief Macro INSTANTIATE_TYPE_INDICES is used to to instantiate a template class with
 * various predefined index configurations.
 *
 * \param CLASS The name of the class template to instantiate.
 * \param T The scalar type (e.g., `double` or `float`) used in the instantiation.
 *
 * \details
 * The macro internally uses `INSTANTIATE_CLASS` to instantiate the provided class template
 * with combinations of `BlackOilOnePhaseIndices`, `BlackOilTwoPhaseIndices`, and
 * `BlackOilVariableAndEquationIndices` (three phase variants). These indices represent
 * different configurations for different physical systems to be simulated.
 *
 * \usage
 * Example usage in a source file:
 * \code
 * #include <opm/simulators/utils/InstantiationIndicesMacros.hpp>
 *
 * INSTANTIATE_TYPE_INDICES(MyClass, double)
 *
 * #if FLOW_INSTANTIATE_FLOAT
 * INSTANTIATE_TYPE_INDICES(MyClass, float)
 * #endif
 * \endcode
 *
 * This will instantiate `MyClass` with all predefined index configurations for both
 * `double` and `float` scalar types (if `FLOW_INSTANTIATE_FLOAT` is defined).
 */

    template<class Scalar>
    using FS = Opm::BlackOilFluidSystem<Scalar, Opm::BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE_CLASS(CLASS, T, ...) \
    template class CLASS<FS<T>,__VA_ARGS__>;

#define INSTANTIATE_TYPE_INDICES(CLASS, T)                                                  \
    INSTANTIATE_CLASS(CLASS,T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)  \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)  \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)  \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,2u>) \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,true,false,0u,0u>)             \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,0u,0u,false,true,0u,0u>)             \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,1u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,1u,0u,false,false,0u,0u>)            \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,false,0u,0u>)            \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<0u,0u,0u,1u,false,true,0u,0u>)             \
    INSTANTIATE_CLASS(CLASS,T,BlackOilVariableAndEquationIndices<1u,0u,0u,0u,true,false,0u,0u>)

#endif
