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

#include "config.h"

#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>

namespace Opm {


VFPProperties::VFPProperties() {
}

VFPProperties::VFPProperties(const VFPInjTable* inj_table, const VFPProdTable* prod_table) {
    if (inj_table != NULL) {
        //FIXME: Implement VFPInjProperties
        OPM_THROW(std::logic_error, "VFPInjProperties not implemented yet");
    }
    if (prod_table != NULL) {
        m_prod.reset(new VFPProdProperties(prod_table));
    }
}

VFPProperties::VFPProperties(const std::map<int, VFPInjTable>& inj_tables,
                             const std::map<int, VFPProdTable>& prod_tables) {
    //FIXME: Implement VFPInjProperties
    OPM_THROW(std::logic_error, "VFPInjProperties not implemented yet");
    m_prod.reset(new VFPProdProperties(prod_tables));
}


VFPProperties::VFPProperties(const std::map<int, VFPInjTable>& inj_tables) {
    //FIXME: Implement VFPInjProperties
    OPM_THROW(std::logic_error, "VFPInjProperties not implemented yet");
}


VFPProperties::VFPProperties(const std::map<int, VFPProdTable>& prod_tables) {
    m_prod.reset(new VFPProdProperties(prod_tables));
}



} //Namespace
