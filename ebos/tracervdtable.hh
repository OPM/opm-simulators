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
 *
 * \copydoc Ewoms::TracerVdTable
 */
#ifndef EWOMS_TRACER_VD_TABLE_HH
#define EWOMS_TRACER_VD_TABLE_HH

#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>

namespace Ewoms {

/*!
 * \brief A class that contains tracer concentration vs depth table
 */
class TracerVdTable : public Opm::SimpleTable
{
public:
    TracerVdTable(const Opm::DeckItem& item)
    {
        this->m_schema.addColumn(Opm::ColumnSchema("DEPTH", Opm::Table::STRICTLY_INCREASING, Opm::Table::DEFAULT_NONE));
        this->m_schema.addColumn(Opm::ColumnSchema("TRACER_CONCENTRATION", Opm::Table::RANDOM, Opm::Table::DEFAULT_NONE));

        Opm::SimpleTable::init(item);
    }

    /*!
     * \brief Return the depth column
     */
    const Opm::TableColumn& getDepthColumn() const
    { return Opm::SimpleTable::getColumn(0); }

    /*!
     * \brief Return the tracer concentration column
     */
    const Opm::TableColumn& getTracerConcentration() const
    { return Opm::SimpleTable::getColumn(1); }
};
}


#endif // TRACERVDTABLE_HH

