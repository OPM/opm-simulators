/*
  Copyright (C) 2013 by Andreas Lauser

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
*/
/*!
 * \file
 * \brief This file contains helper classes for the material laws.
 *
 * These classes specify the information which connects fluid systems
 * and the fluid-matrix interaction laws. This includes stuff like the
 * index of the wetting and non-wetting phase, etc.
 */
#ifndef OPM_MATERIAL_TRAITS_HPP
#define OPM_MATERIAL_TRAITS_HPP

namespace Opm {
/*!
 * \ingroup material
 *
 * \brief A generic traits class which does not provide any indices.
 *
 * This traits class is intended to be used by the NullMaterial
 */
template <class ScalarT, int numPhasesV>
class NullMaterialTraits
{
public:
    //! The type used for scalar floating point values
    typedef ScalarT Scalar;

    //! The number of fluid phases
    static const int numPhases = numPhasesV;
};

/*!
 * \ingroup material
 *
 * \brief A generic traits class for two-phase material laws.
 */
template <class ScalarT, int wettinPhaseIdxV, int nonWettingasPhaseIdxV>
class TwoPhaseMaterialTraits
{
public:
    //! The type used for scalar floating point values
    typedef ScalarT Scalar;

    //! The number of fluid phases
    static const int numPhases = 2;

    //! The index of the wetting phase
    static const int  wettingPhaseIdx = wettinPhaseIdxV;

    //! The index of the non-wetting phase
    static const int nonWettingPhaseIdx = nonWettingasPhaseIdxV;

    // some safety checks...
    static_assert(wettingPhaseIdx != nonWettingPhaseIdx,
                  "wettingPhaseIdx and nonWettingPhaseIdx must be different");
};

/*!
 * \ingroup material
 *
 * \brief A generic traits class for three-phase material laws.
 */
template <class ScalarT, int wettingPhaseIdxV, int nonWettingasPhaseIdxV, int gasPhaseIdxV>
class ThreePhaseMaterialTraits
{
public:
    //! The type used for scalar floating point values
    typedef ScalarT Scalar;

    //! The number of fluid phases
    static const int numPhases = 3;

    //! The index of the wetting liquid phase
    static const int wettingPhaseIdx = wettingPhaseIdxV;

    //! The index of the non-wetting liquid phase
    static const int nonWettingPhaseIdx = nonWettingasPhaseIdxV;

    //! The index of the gas phase (i.e., the least wetting phase)
    static const int gasPhaseIdx = gasPhaseIdxV;

    // some safety checks...
    static_assert(0 <= wettingPhaseIdx && wettingPhaseIdx < numPhases,
                  "wettingPhaseIdx is out of range");
    static_assert(0 <= nonWettingPhaseIdx && nonWettingPhaseIdx < numPhases,
                  "nonWettingPhaseIdx is out of range");
    static_assert(0 <= gasPhaseIdx && gasPhaseIdx < numPhases,
                  "gasPhaseIdx is out of range");

    static_assert(wettingPhaseIdx != nonWettingPhaseIdx,
                  "wettingPhaseIdx and nonWettingPhaseIdx must be different");
    static_assert(wettingPhaseIdx != gasPhaseIdx,
                  "wettingPhaseIdx and gasPhaseIdx must be different");
    static_assert(nonWettingPhaseIdx != gasPhaseIdx,
                  "nonWettingPhaseIdx and gasPhaseIdx must be different");
};
} // namespace Opm

#endif
