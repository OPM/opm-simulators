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
    class WellProdIndexCalculator
    {
    public:
        explicit WellProdIndexCalculator(const Well& well);

        double connectionProdIndStandard(const std::size_t connIdx,
                                         const double      connMobility) const;

        std::size_t numConnections() const
        {
            return this->standardConnFactors_.size();
        }

    private:
        std::vector<double> standardConnFactors_{};
    };

    std::vector<double>
    connectionProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                              const std::vector<double>&     connMobility);

    double wellProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                               const std::vector<double>&     connMobility);
} // namespace Opm

#endif // OPM_WELLPRODINDEXCALCULATOR_HEADER_INCLUDED
