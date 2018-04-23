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

#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>

#include <map>

namespace Opm {

class VFPProdProperties;
class VFPInjProperties;

/**
 * A thin wrapper class that holds one VFPProdProperties and one
 * VFPInjProperties object.
 */
class VFPProperties {
public:
    VFPProperties();

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_table  A *single* VFPINJ table or NULL (no table)
     * @param prod_table A *single* VFPPROD table or NULL (no table)
     */
    explicit VFPProperties(const VFPInjTable* inj_table, const VFPProdTable* prod_table);

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     * @param prod_tables A map of different VFPPROD tables.
     */
    VFPProperties(const std::map<int, std::shared_ptr<const VFPInjTable> >& inj_tables,
                  const std::map<int, std::shared_ptr<const VFPProdTable> >& prod_tables);

    /**
     * Returns the VFP properties for injection wells
     */
    const VFPInjProperties* getInj() const {
        return m_inj.get();
    }

    /**
     * Returns the VFP properties for production wells
     */
    const VFPProdProperties* getProd() const {
        return m_prod.get();
    }

private:
    std::shared_ptr<VFPInjProperties> m_inj;
    std::shared_ptr<VFPProdProperties> m_prod;
};


} //Namespace

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
