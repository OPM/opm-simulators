/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_WELLPRODINDEXCALCULATOR_HEADER_INCLUDED
#define OPM_WELLPRODINDEXCALCULATOR_HEADER_INCLUDED

#include <cstddef>
#include <vector>

namespace Opm {
    class Well;
} // namespace Opm

namespace Opm {

    /// Collect per-connection static information to enable calculating
    /// connection-level or well-level productivity index values when
    /// incorporating dynamic phase mobilities.
    class WellProdIndexCalculator
    {
    public:
        /// Constructor
        ///
        /// \param[in] well Individual well for which to collect
        ///   per-connection static data.
        explicit WellProdIndexCalculator(const Well& well);

        /// Reinitialization operation
        ///
        /// Needed to repopulate the internal data members in case of
        /// changes to the Well's properties, e.g., as a result of the
        /// Well's CTFs being rescaled due to WELPI.
        ///
        /// \param[in] well Individual well for which to collect
        ///   per-connection static data.
        void reInit(const Well& well);

        /// Compute connection-level steady-state productivity index value
        /// using dynamic phase mobility.
        ///
        /// \param[in] connIdx Linear connection index.  Must be in the
        ///   range 0..numConnections() - 1.
        ///
        /// \param[in] connMobility Phase mobility at connection \p connIdx.
        ///   Typically derived from dynamic flow state conditions in cell
        ///   intersected by well's connection \p connIdx.
        ///
        /// \return Connection-level steady-state productivity index.
        double connectionProdIndStandard(const std::size_t connIdx,
                                         const double      connMobility) const;

        /// Number of connections in this well.
        ///
        /// Used primarily for consistency checks.
        std::size_t numConnections() const
        {
            return this->standardConnFactors_.size();
        }

    private:
        /// Static, per-connection multiplicative PI factors.
        ///
        /// Corresponds to the well's connection transmissibility factors,
        /// multiplied by a ratio of logarithms if the well has an explicit,
        /// positive drainage radius.
        std::vector<double> standardConnFactors_{};
    };

    /// Compute connection-level productivity index values for all
    /// connections in a well.
    ///
    /// \param[in] wellPICalc Productivity index calculator.
    ///
    /// \param[in] connMobility Phase mobility for each connection.
    ///   Typically derived from dynamic flow state conditions in cells
    ///   intersected by well's connections.  Must have one value for each
    ///   \code wellPICalc.numConnections() \endcode well connection.
    ///
    /// \return Connection-level steady-state productivity index values for
    ///   all connections.
    std::vector<double>
    connectionProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                              const std::vector<double>&     connMobility);

    /// Compute well-level productivity index value.
    ///
    /// \param[in] wellPICalc Productivity index calculator.
    ///
    /// \param[in] connMobility Phase mobility for each connection.
    ///   Typically derived from dynamic flow state conditions in cells
    ///   intersected by well's connections.  Must have one value for each
    ///   \code wellPICalc.numConnections() \endcode well connection.
    ///
    /// \return Well-level steady-state productivity index value.
    double wellProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                               const std::vector<double>&     connMobility);
} // namespace Opm

#endif // OPM_WELLPRODINDEXCALCULATOR_HEADER_INCLUDED
