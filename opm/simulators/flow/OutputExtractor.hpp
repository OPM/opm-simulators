// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/*!
 * \file
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_OUTPUT_EXTRACTORS_HPP
#define OPM_OUTPUT_EXTRACTORS_HPP

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <array>
#include <variant>
#include <vector>

namespace Opm::detail {

//! \brief Wrapping struct holding types used for element-level data extraction.
template<class TypeTag>
struct Extractor
{
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidState = typename IntensiveQuantities::FluidState;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numPhases = FluidSystem::numPhases;

    //! \brief Struct holding hysteresis parameters.
    struct HysteresisParams
    {
        Scalar somax{}; //!< Max oil saturation
        Scalar swmax{}; //!< Max water saturation
        Scalar swmin{}; //!< Min water saturation
        Scalar sgmax{}; //!< Max gas saturation
        Scalar shmax{}; //!< Max something
        Scalar somin{}; //!< Min oil saturation
    };

    //! \brief Context passed to extractor functions.
    struct Context
    {
        unsigned globalDofIdx; //!< Global degree-of-freedom index
        unsigned pvtRegionIdx; //!< pvt region for dof
        int episodeIndex;      //!< Current report step
        const FluidState& fs;  //!< Fluid state for cell
        const IntensiveQuantities& intQuants; //!< Intensive quantities for cell
        const HysteresisParams& hParams; //!< Hysteresis parameters for cell
    };

    /// Callback for extractors handling their own assignements
    using AssignFunc = std::function<void(const Context&)>;

    /// Callback for extractors assigned to a scalar buffer
    /// Return value to store in buffer
    using ScalarFunc = std::function<Scalar(const Context&)>;

    /// Callback for extractors assigned to a phase buffer
    /// Returns value to store in buffer for requested phase
    using PhaseFunc = std::function<Scalar(const unsigned /*phase*/, const Context&)>;

    using ScalarBuffer = std::vector<Scalar>; //!< A scalar buffer
    using PhaseArray = std::array<ScalarBuffer,numPhases>; //!< An array of buffers, one for each phase

    //! \brief A scalar extractor descriptor.
    struct ScalarEntry
    {
        ScalarBuffer* data; //!< Buffer to store data in
        ScalarFunc extract; //!< Function to call for extraction
    };

    //! \brief A phase buffer extractor descriptor.
    struct PhaseEntry
    {
        PhaseArray* data; //!< Array of buffers to store data in
        PhaseFunc extract; //!< Function to call for extraction
    };

    //! \brief Descriptor for extractors
    struct Entry
    {
        std::variant<AssignFunc, ScalarEntry, PhaseEntry> data; //!< Extractor
        bool condition = true; //!< Additional condition for enabling extractor
    };
};

} // namespace Opm::detail

#endif // OPM_OUTPUT_EXTRACTORS_HPP
