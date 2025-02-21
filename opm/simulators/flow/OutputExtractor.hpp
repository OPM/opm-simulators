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

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/Visitor.hpp>

#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <set>
#include <variant>
#include <vector>

#include <fmt/format.h>

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

    //! \brief Obtain vector of active extractors from an array of extractors.
    template<std::size_t size>
    static std::vector<Entry> removeInactive(std::array<Entry,size>& input)
    {
        // Setup active extractors
        std::vector<Entry> filtered_extractors;
        filtered_extractors.reserve(input.size());
        std::copy_if(std::move_iterator(input.begin()),
                     std::move_iterator(input.end()),
                     std::back_inserter(filtered_extractors),
                     [](const Entry& e)
                     {
                         if (!e.condition) {
                            return false;
                         }
                         return std::visit(VisitorOverloadSet{
                                               [](const AssignFunc&)
                                               {
                                                   return true;
                                               },
                                               [](const ScalarEntry& v)
                                               {
                                                   return !v.data->empty();
                                               },
                                               [](const PhaseEntry& v)
                                               {
                                                   return std::any_of(v.data->begin(),
                                                                      v.data->end(),
                                                                      [](const auto& ve)
                                                                      { return !ve.empty(); });
                                               }
                                           }, e.data);
                     });

        return filtered_extractors;
    }

    //! \brief Process the given extractor entries
    //! \param ectx Context for extractors
    //! \param extractors List of extractors to process
    static void process(const Context& ectx,
                        const std::vector<Entry>& extractors)
    {
        std::for_each(extractors.begin(), extractors.end(),
                      [&ectx](const auto& entry)
                      {
                          std::visit(VisitorOverloadSet{
                              [&ectx](const ScalarEntry& v)
                              {
                                  auto& array = *v.data;
                                  array[ectx.globalDofIdx] = v.extract(ectx);
                                  Valgrind::CheckDefined(array[ectx.globalDofIdx]);
                              },
                              [&ectx](const PhaseEntry& v)
                              {
                                  std::for_each(v.data->begin(), v.data->end(),
                                                [phaseIdx = 0, &ectx, &v](auto& array) mutable
                                                {
                                                    if (!array.empty()) {
                                                        array[ectx.globalDofIdx] = v.extract(phaseIdx, ectx);
                                                        Valgrind::CheckDefined(array[ectx.globalDofIdx]);
                                                    }
                                                    ++phaseIdx;
                                                });
                              },
                              [&ectx](const typename Extractor::AssignFunc& extract)
                              { extract(ectx); },
                          }, entry.data);
                      });
    }
};

//! \brief Wrapping struct holding types used for block-level data extraction.
template<class TypeTag>
struct BlockExtractor
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidState = typename IntensiveQuantities::FluidState;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;

    //! \brief Context passed to element extractor functions.
    struct Context
    {
        unsigned globalDofIdx; //!< Global degree-of-freedom index
        unsigned dofIdx;
        const FluidState& fs;  //!< Fluid state for cell
        const IntensiveQuantities& intQuants; //!< Intensive quantities for cell
        const ElementContext& elemCtx;
    };

    /// Callback for extractors handling their own assignements
    using AssignFunc = std::function<void(const Context&)>;

    /// Callback for extractors assigned to a scalar buffer
    /// Return value to store in buffer
    using ScalarFunc = std::function<Scalar(const Context&)>;

    /// Callback for extractors assigned to a phase buffer
    /// Returns value to store in buffer for requested phase
    using PhaseFunc = std::function<Scalar(const unsigned /*phase*/, const Context&)>;

    struct ScalarEntry
    {
        /// A single name or a list of names for the keyword
        std::variant<std::string_view, std::vector<std::string_view>> kw;

        /// Associated extraction lamda
        ScalarFunc extract;
    };

    struct PhaseEntry
    {
        /// One or two lists of names for the keyword for each phase
        std::variant<std::array<std::string_view, numPhases>,
                     std::array<std::array<std::string_view, numPhases>, 2>> kw;

        /// Associated extraction lambda
        PhaseFunc extract;
    };

    //! \brief Descriptor for extractors
    using Entry = std::variant<ScalarEntry, PhaseEntry>;

    //! \brief Descriptor for extractor execution.
    struct Exec
    {
        //! \brief Move constructor.
        Exec(double* d, ScalarFunc&& e)
            : data(d), extract(std::move(e))
        {}

        double* data; //!< Where to store output data
        ScalarFunc extract; //!< Extraction function to call
    };

    /// A map of extraction executors, keyed by cartesian cell index
    using ExecMap = std::unordered_map<int, std::vector<Exec>>;

    //! \brief Setup an extractor executor map from a map of evaluations to perform.
    template<std::size_t size>
    static ExecMap setupExecMap(std::map<std::pair<std::string, int>, double>& blockData,
                                const std::array<Entry,size>& handlers)
    {
        using PhaseViewArray = std::array<std::string_view, numPhases>;
        using StringViewVec = std::vector<std::string_view>;

        ExecMap extractors;

        std::for_each(
            blockData.begin(),
            blockData.end(),
            [&handlers, &extractors, &blockData](auto& bd_info)
            {
                unsigned phase{};
                const auto& [key, cell] = bd_info.first;
                const auto& handler_info =
                    std::find_if(
                        handlers.begin(),
                        handlers.end(),
                        [&kw_name = bd_info.first.first, &phase](const auto& handler)
                        {
                           // Extract list of keyword names from handler
                           const auto gen_handlers =
                              std::visit(VisitorOverloadSet{
                                            [](const ScalarEntry& entry)
                                            {
                                                return std::visit(VisitorOverloadSet{
                                                                      [](const std::string_view& kw) -> StringViewVec
                                                                      {
                                                                          return {kw};
                                                                      },
                                                                      [](const StringViewVec& kws) -> StringViewVec
                                                                      { return kws; }
                                                                  }, entry.kw);
                                            },
                                            [](const PhaseEntry& entry)
                                            {
                                                return std::visit(VisitorOverloadSet{
                                                                      [](const PhaseViewArray& data)
                                                                      {
                                                                          return StringViewVec{data.begin(), data.end()};
                                                                      },
                                                                      [](const std::array<PhaseViewArray,2>& data)
                                                                      {
                                                                          StringViewVec res;
                                                                          res.reserve(2*numPhases);
                                                                          res.insert(res.end(),
                                                                                     data[0].begin(),
                                                                                     data[0].end());
                                                                          res.insert(res.end(),
                                                                                     data[1].begin(),
                                                                                     data[1].end());
                                                                          return res;
                                                                      }
                                                                  }, entry.kw);
                                            }
                                        }, handler);

                            const auto found_handler =
                                std::find(gen_handlers.begin(), gen_handlers.end(), kw_name);
                            if (found_handler != gen_handlers.end()) {
                                phase = std::distance(gen_handlers.begin(), found_handler) % numPhases;
                            }
                            return found_handler != gen_handlers.end();
                        }
                    );

                if (handler_info != handlers.end()) {
                    extractors[cell - 1].emplace_back(
                        &bd_info.second,
                        std::visit(VisitorOverloadSet{
                                       [](const ScalarEntry& e)
                                       {
                                           return e.extract;
                                       },
                                       [phase](const PhaseEntry& e) -> ScalarFunc
                                       {
                                           return [phase, extract = e.extract]
                                                  (const Context& ectx)
                                                  {
                                                      static constexpr auto phaseMap = std::array{
                                                          waterPhaseIdx,
                                                          oilPhaseIdx,
                                                          gasPhaseIdx,
                                                      };
                                                      return extract(phaseMap[phase], ectx);
                                                  };
                                       }
                                   }, *handler_info)
                     );
                }
                else {
                    OpmLog::warning("Unhandled output keyword",
                                    fmt::format("Keyword '{}' is unhandled for output "
                                                "to summary file.", key));
                }
            }
        );

        return extractors;
    }

    //! \brief Process a list of block extractors.
    static void process(const std::vector<Exec>& blockExtractors,
                        const Context& ectx)
    {
        std::for_each(blockExtractors.begin(), blockExtractors.end(),
                      [&ectx](auto& bdata)
                      { *bdata.data = bdata.extract(ectx); });
    }
};

} // namespace Opm::detail

#endif // OPM_OUTPUT_EXTRACTORS_HPP
