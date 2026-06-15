/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_FLOWEXP_GRAPH_WELL_MODEL_HEADER_INCLUDED
#define OPM_FLOWEXP_GRAPH_WELL_MODEL_HEADER_INCLUDED

#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <flowexperimental/graphwell/GraphMultisegmentWell.hpp>

namespace Opm {

/// \brief Well model that substitutes the GraphWell multisegment well for the
/// production MultisegmentWell. Standard wells are left untouched.
template<typename TypeTag>
class FlowExpGraphWellModel : public BlackoilWellModel<TypeTag>
{
public:
    using Base = BlackoilWellModel<TypeTag>;
    using WellInterfacePtr = typename Base::WellInterfacePtr;

    using Base::Base;   // inherit constructors

    WellInterfacePtr
    createWellPointer(const int wellID, const int report_step) const override
    {
        const auto is_multiseg = this->wells_ecl_[wellID].isMultiSegment();
        if (this->param_.use_multisegment_well_ && is_multiseg) {
            return this->template createTypedWellPointer<GraphMultisegmentWell<TypeTag>>(wellID, report_step);
        }
        return Base::createWellPointer(wellID, report_step);
    }
};

} // namespace Opm

#endif // OPM_FLOWEXP_GRAPH_WELL_MODEL_HEADER_INCLUDED
