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
 * \brief This file contains the default flux module of the blackoil model.
 *
 * It is neccessary to accomodate the extensions of the black-oil model.
 */
#ifndef EWOMS_BLACK_OIL_DARCY_FLUX_MODULE_HH
#define EWOMS_BLACK_OIL_DARCY_FLUX_MODULE_HH

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/common/darcyfluxmodule.hh>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

template <class TypeTag>
class BlackOilDarcyExtensiveQuantities;

/*!
 * \ingroup FluxModules
 * \brief Provides a Darcy flux module for the blackoil model
 */
template <class TypeTag>
struct BlackOilDarcyFluxModule
{
    typedef DarcyIntensiveQuantities<TypeTag> FluxIntensiveQuantities;
    typedef BlackOilDarcyExtensiveQuantities<TypeTag> FluxExtensiveQuantities;
    typedef DarcyBaseProblem<TypeTag> FluxBaseProblem;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};


/*!
 * \ingroup FluxModules
 * \brief Specifies the extensive quantities for the black-oil model if using Darcy relation.
 *
 * This class basically forwards everything to the default Darcy flux module and adds a
 * few methods needed by the extensions of the black-oil model. (i.e. the solvent and the
 * polymer extensions.)
 */
template <class TypeTag>
class BlackOilDarcyExtensiveQuantities : public DarcyExtensiveQuantities<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    /*!
     * \brief Update the extensive quantities which are specific to the solvent extension
     * of the black-oil model.
     */
    void updateSolvent(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        asImp_().updateVolumeFluxPerm(elemCtx,
                                      scvfIdx,
                                      timeIdx);

    }

    void updatePolymer(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    { asImp_().updateShearMultipliersPerm(elemCtx, scvfIdx, timeIdx); }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
};

} // namespace Opm

#endif
