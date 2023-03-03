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
 * \copydoc Opm::EclBaseAquiferModel
 */
#ifndef EWOMS_ECL_BASE_AQUIFER_MODEL_HH
#define EWOMS_ECL_BASE_AQUIFER_MODEL_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/output/data/Aquifer.hpp>

#include <exception>
#include <stdexcept>
#include <vector>

namespace Opm {

/*!
 * \ingroup EclBaseAquiferModel
 *
 * \brief The base class which specifies the API of aquifer models.
 *
 * This class only provides the API for the actual aquifer model, it does not do
 * anything on its own.
 */
template <class TypeTag>
class EclBaseAquiferModel
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;

public:
    EclBaseAquiferModel(Simulator& simulator)
        : simulator_(simulator)
    {}

    /*!
     * \brief Called once the problem has been fully initialized and the initial
     *        condition has been applied.
     */
    void initialSolutionApplied()
    { }

    /*!
     * \brief Called if aquifers are being initialized from values retrieved
     *        from a restart file.
     *
     * \param[in] aquiferSoln Set of aquifer-related initial values, mostly
     *        pertaining to analytic aquifers.  Contains at minimum the
     *        aquifer pressure and the base run's total produced liquid
     *        volume from the model's aquifers.
     */
    void initFromRestart(const data::Aquifers&)
    {
        throw std::logic_error {
            "Initialization from restart data not supported "
            "for base aquifer model"
        };
    }

    /*!
     * \brief This method is called when a new episode (report step) starts.
     */
    void beginEpisode()
    { }

    /*!
     * \brief This method is called when a new time step (substep) starts.
     */
    void beginTimeStep()
    { }

    /*!
     * \brief This method is called before each Newton-Raphson iteration.
     */
    void beginIteration()
    { }

    /*!
     * \brief Add the water which enters or leaves the reservoir due to aquifiers.
     */
    template <class Context>
    void addToSource(RateVector&,
                     const Context&,
                     unsigned,
                     unsigned) const
    { }

    /*!
     * \brief This method is called after each Newton-Raphson successful iteration.
     *
     * I.e., no exceptions were thrown during the linearization and linear solution
     * procedures.
     */
    void endIteration()
    { }

    /*!
     * \brief This method is called after each successful time step (substep).
     *
     * I.e., all iterations of the time step were successful and the Newton-Raphson
     * algorithm converged.
     */
    void endTimeStep()
    { }

    /*!
     * \brief This method is called once an episode (report step) has been finished
     *        successfully.
     */
    void endEpisode()
    { }

    /*!
     * \brief Write the internal state of the aquifer model to disk using an ad-hoc file
     *        format.
     */
    template <class Restarter>
    void serialize(Restarter&)
    { }

    /*!
     * \brief Load the internal state of the aquifer model to disk using an ad-hoc file
     *        format.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    { }


    data::Aquifers aquiferData() const
    { return data::Aquifers{}; }


protected:
    Simulator& simulator_;
};

} // namespace Opm

#endif
