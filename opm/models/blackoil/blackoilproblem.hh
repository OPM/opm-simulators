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
 * \copydoc Opm::BlackOilProblem
 */
#ifndef EWOMS_BLACKOIL_PROBLEM_HH
#define EWOMS_BLACKOIL_PROBLEM_HH

#include "blackoilproperties.hh"

#include <opm/models/common/multiphasebaseproblem.hh>

#include <opm/material/common/Unused.hpp>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 * \brief Base class for all problems which use the black-oil model.
 */
template<class TypeTag>
class BlackOilProblem : public MultiPhaseBaseProblem<TypeTag>
{
private:
    typedef MultiPhaseBaseProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     *
     * \param simulator The manager object of the simulation
     */
    BlackOilProblem(Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Returns the maximum value of the gas dissolution factor at the current time
     *        for a given degree of freedom.
     *
     * This is required for the DRSDT keyword.
     */
    Scalar maxGasDissolutionFactor(unsigned timeIdx OPM_UNUSED, unsigned globalDofIdx OPM_UNUSED) const
    { return std::numeric_limits<Scalar>::max()/2; }

    /*!
     * \brief Returns the maximum value of the oil vaporization factor at the current
     *        time for a given degree of freedom.
     *
     * This is required for the DRVDT keyword.
     */
    Scalar maxOilVaporizationFactor(unsigned timeIdx OPM_UNUSED, unsigned globalDofIdx OPM_UNUSED) const
    { return std::numeric_limits<Scalar>::max()/2; }

    /*!
     * \brief Returns the maximum value of the oil saturation seen at the current time
     *        for a given degree of freedom.
     *
     * This is required for the VAPPARS keyword.
     */
    Scalar maxOilSaturation(unsigned globalDofIdx OPM_UNUSED) const
    { return 1.0; }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    { return 0; }

    /*!
     * \brief Returns the index of the relevant region for saturation functions
     */
    template <class Context>
    unsigned satnumRegionIndex(const Context& context OPM_UNUSED,
                               unsigned spaceIdx OPM_UNUSED,
                               unsigned timeIdx OPM_UNUSED) const
    { return 0; }

    /*!
     * \brief Returns the index of the relevant region for solvent mixing functions
     */
    template <class Context>
    unsigned miscnumRegionIndex(const Context& context OPM_UNUSED,
                                unsigned spaceIdx OPM_UNUSED,
                                unsigned timeIdx OPM_UNUSED) const
    { return 0; }

    /*!
     * \brief Returns the index of the relevant region for polymer mixing functions
     */
    template <class Context>
    unsigned plmixnumRegionIndex(const Context& context OPM_UNUSED,
                                 unsigned spaceIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED) const
    { return 0; }

    /*!
     * \brief Returns the compressibility of the porous medium of a cell
     */
    template <class Context>
    Scalar rockCompressibility(const Context& context OPM_UNUSED,
                               unsigned spaceIdx OPM_UNUSED,
                               unsigned timeIdx OPM_UNUSED) const
    { return 0.0; }

    /*!
     * \brief Returns the reference pressure for rock the compressibility of a cell
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context OPM_UNUSED,
                                 unsigned spaceIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED) const
    { return 1e5; }

    /*!
     * \brief Returns the reference temperature
     *
     * This is only relevant for temperature dependent quantities, in particular those
     * needed by the module for energy conservation.
     */
    Scalar referenceTemperature() const
    { return 273.15 + 15.56; /* [K] */ }

    /*!
     * \brief Returns the porosity multiplier due to water-induced rock compaction
     *
     * This is a somewhat exotic feature. Most likely you will not need to touch this
     * method.
     */
    template <class Evaluation>
    Scalar rockCompPoroMultiplier(const IntensiveQuantities& intQuants OPM_UNUSED,
                                  unsigned globalSpaceIdx OPM_UNUSED) const
    { return 1.0; }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Opm

#endif
