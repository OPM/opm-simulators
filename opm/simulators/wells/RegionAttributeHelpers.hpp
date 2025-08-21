/*
  Copyright 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2017, IRIS
  Copyright 2017, Equinor

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_REGIONATTRIBUTEHELPERS_HPP_HEADER_INCLUDED
#define OPM_REGIONATTRIBUTEHELPERS_HPP_HEADER_INCLUDED

#include <opm/grid/utility/RegionMapping.hpp>

#include <dune/grid/common/gridenums.hh>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm {
    namespace RegionAttributeHelpers {
        /**
         * Convenience tools for processing region
         * spesific attributes
         */
            namespace Select {
                template <class RegionID, bool>
                struct RegionIDParameter
                {
                    using type =
                        typename std::remove_reference<RegionID>::type &;
                };

                template <class RegionID>
                struct RegionIDParameter<RegionID, true>
                {
                    using type = RegionID;
                };
            } // Select

            /**
             * \brief Computes the temperature, pressure, and counter increment.
             *
             * In a parallel run only cells owned contribute to the cell average.
             * \tparam is_parallel Whether this is a parallel run.
             */
            template<class Scalar, bool is_parallel>
            struct AverageIncrementCalculator
            {
                /**
                 * \brief Computes the temperature, pressure, and counter increment.
                 * \param pressure    The pressure.
                 * \param temperature The temperature.
                 * \param rs          The rs.
                 * \param rv          The rv.
                 * \param cell        The current cell index.
                 * \param ownership   A vector indicating whether a cell is owned
                 *                    by this process (value 1), or not (value 0).
                 * \param cell        The cell index.
                 */
                std::tuple<Scalar, Scalar, Scalar, Scalar, int>
                operator()(const std::vector<Scalar>& pressure,
                           const std::vector<Scalar>& temperature,
                           const std::vector<Scalar>& rs,
                           const std::vector<Scalar>& rv,
                           const std::vector<Scalar>& ownership,
                           std::size_t cell){
                    if ( ownership[cell] )
                    {
                        return std::make_tuple(pressure[cell],
                                               temperature[cell],
                                               rs[cell],
                                               rv[cell],
                                               1);
                    }
                    else
                    {
                        return std::make_tuple(0, 0, 0, 0, 0);
                    }
                }
            };
            template<class Scalar>
            struct AverageIncrementCalculator<Scalar, false>
            {
                std::tuple<Scalar, Scalar, Scalar, Scalar, int>
                operator()(const std::vector<Scalar>& pressure,
                           const std::vector<Scalar>& temperature,
                           const std::vector<Scalar>& rs,
                           const std::vector<Scalar>& rv,
                           const std::vector<Scalar>&,
                           std::size_t cell){
                    return std::make_tuple(pressure[cell],
                                           temperature[cell],
                                           rs[cell],
                                           rv[cell],
                                           1);
                }
            };
            /**
             * Provide mapping from Region IDs to user-specified collection
             * of per-region attributes.
             *
             * \tparam RegionId Region identifier type.  Must be hashable by
             * \code std::hash<> \endcode.  Typically a built-in integer
             * type--e.g., \c int.
             *
             * \tparam Attributes User-defined type that represents
             * collection of attributes that have meaning in a per-region
             * aggregate sense.  Must be copy-constructible.
             */
            template <typename RegionId, class Attributes>
            class RegionAttributes
            {
            public:
                /**
                 * Expose \c RegionId as a vocabulary type for use in query
                 * methods.
                 */
                using RegionID =
                    typename Select::RegionIDParameter
                    <RegionId, std::is_integral<RegionId>::value>::type;

                using ID =
                    typename std::remove_reference<RegionId>::type;

                /**
                 * Aggregate per-region attributes along with region's
                 * representative cell.
                 */
                struct Value {
                    explicit Value(const Attributes& attr)
                        : attr_(attr)
                        , cell_(-1)
                    {}

                    Attributes attr_;
                    int        cell_;
                };

                using AttributeMap =
                    std::unordered_map<ID, std::unique_ptr<Value>>;


                /**
                 * Constructor.
                 *
                 * \tparam RMap Class type that implements the RegionMapping
                 * protocol.  Typically an instantiation of \code
                 * Opm::RegionMapping<> \endcode.
                 *
                 * \param[in] rmap Specific region mapping that provides
                 * reverse lookup from regions to cells.
                 *
                 * \param[in] attr Pre-constructed initialiser for \c
                 * Attributes.
                 */
                template <class RMap>
                RegionAttributes(const RMap&       rmap,
                                 const Attributes& attr)
                {
                    using VT = typename AttributeMap::value_type;

                    for (const auto& r : rmap.activeRegions()) {
                        auto v = std::make_unique<Value>(attr);

                        const auto stat = attr_.insert(VT(r, std::move(v)));

                        if (stat.second) {
                            // New value inserted.
                            const auto& cells = rmap.cells(r);

                            assert (! cells.empty());

                            // Region's representative cell.
                            stat.first->second->cell_ = cells[0];
                        }
                    }
                }

                /**
                 * Retrieve representative cell in region.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Representative cell in region \p reg.
                 */
                int cell(const RegionID reg) const
                {
                    return this->find(reg).cell_;
                }

                bool has(const RegionID reg) const
                {
                    return this->attr_.find(reg) != this->attr_.end();
                }

                void insert(const RegionID r, const Attributes& attr)
                {
                    auto [pos, inserted] = this->attr_.try_emplace(r, std::make_unique<Value>(attr));
                    if (inserted) {
                        pos->second->cell_ = -1; // NOT -1.0 -- "cell_" is 'int'
                    }
                }

                /**
                 * Request read-only access to region's attributes.
                 *
                *
                 * \return Read-only access to all regions attributes.
                 */
                const AttributeMap& attributes() const
                {
                    return attr_;
                }


                /**
                 * Request read-only access to region's attributes.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Read-only access to region \p reg's per-region
                 * attributes.
                 */
                const Attributes& attributes(const RegionID reg) const
                {
                    return this->find(reg).attr_;
                }

                /**
                 * Request modifiable access to region's attributes.
                 *
                 * \param[in] reg Specific region.
                 *
                 * \return Read-write access to region \p reg's per-region
                 * attributes.
                 */
                Attributes& attributes(const RegionID reg)
                {
                    return this->find(reg).attr_;
                }

            private:

                AttributeMap attr_;

                /**
                 * Read-only access to region's properties.
                 */
                const Value& find(const RegionID reg) const
                {
                    const auto& i = attr_.find(reg);

                    if (i == attr_.end()) {
                        throw std::invalid_argument("Unknown region ID");
                    }

                    return *i->second;
                }

                /**
                 * Read-write access to region's properties.
                 */
                Value& find(const RegionID reg)
                {
                    const auto& i = attr_.find(reg);

                    if (i == attr_.end()) {
                        throw std::invalid_argument("Unknown region ID");
                    }

                    return *i->second;
                }
            };
    } // namespace RegionAttributesHelpers
} // namespace Opm

#endif  /* OPM_REGIONATTRIBUTEHELPERS_HPP_HEADER_INCLUDED */
