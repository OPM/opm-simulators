/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFF_VFPPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPPROPERTIES_HPP_

#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>

#include <map>

namespace Opm {

class VFPInjTable;
class VFPProdTable;

/**
 * A thin wrapper class that holds one VFPProdProperties and one
 * VFPInjProperties object.
 */
class VFPProperties {
public:
    VFPProperties() = default;


    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     * @param prod_tables A map of different VFPPROD tables.
     */

    VFPProperties(const std::vector<std::reference_wrapper<const VFPInjTable>>& inj_tables,
                  const std::vector<std::reference_wrapper<const VFPProdTable>>& prod_tables)
    {
        for (const auto& vfpinj : inj_tables)
            this->m_inj.addTable( vfpinj );

        for (const auto& vfpprod : prod_tables)
            this->m_prod.addTable( vfpprod );
    };

    /**
     * Returns the VFP properties for injection wells
     */
    const VFPInjProperties* getInj() const {
        return &m_inj;
    }

    /**
     * Returns the VFP properties for production wells
     */
    const VFPProdProperties* getProd() const {
        return &m_prod;
    }

private:
    VFPInjProperties m_inj;
    VFPProdProperties m_prod;
};


} //Namespace

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
