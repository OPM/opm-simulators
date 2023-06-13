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
 * \copydoc Opm::EclThresholdPressure
 */
#ifndef EWOMS_ECL_GENERIC_THRESHOLD_PRESSURE_HH
#define EWOMS_ECL_GENERIC_THRESHOLD_PRESSURE_HH

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/LookUpData.hh>


#include <vector>

namespace Opm {

class EclipseState;
template<typename Grid, typename GridView> class LookUpData;
template<typename Grid, typename GridView> class LookUpCartesianData;

template<class Grid, class GridView, class ElementMapper, class Scalar>
class EclGenericThresholdPressure {
public:
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using LookUpData = Opm::LookUpData<Grid,GridView>;
    using LookUpCartesianData = Opm::LookUpCartesianData<Grid,GridView>;

    EclGenericThresholdPressure(const CartesianIndexMapper& cartMapper,
                                const GridView& gridView,
                                const ElementMapper& elementMapper,
                                const EclipseState& eclState);

    /*!
     * \brief Returns the theshold pressure [Pa] for the intersection between two elements.
     *
     * This is tailor made for the E100 threshold pressure mechanism and it is thus quite
     * a hack: First of all threshold pressures in general are unphysical, and second,
     * they should be different for the fluid phase but are not. Anyway, this seems to be
     * E100's way of doing things, so we do it the same way.
     */
    Scalar thresholdPressure(int elem1Idx, int elem2Idx) const;

    /*!
     * \brief Return the raw array with the threshold pressures
     *
     * This is used for the restart capability.
     */
    const std::vector<Scalar>& data() const
    { return thpres_; }

    /*!
     * \brief Set the threshold pressures from a raw array
     *
     * This is used for the restart capability.
     */
    void setFromRestart(const std::vector<Scalar>& values)
    { thpres_ = values; }

    //! \brief Returns a fully expanded vector for restart file writing.
    //! \details Returns the union of explicitly configured entries and defaulted values.
    std::vector<Scalar> getRestartVector() const;

protected:
    /*!
     * \brief Actually compute the threshold pressures over a face as a pre-compute step.
     */
    void finishInit();

    // internalize the threshold pressures which where explicitly specified via the
    // THPRES keyword.
    void applyExplicitThresholdPressures_();

    void configureThpresft_();

    void logPressures();

    const CartesianIndexMapper& cartMapper_;
    const GridView& gridView_;
    const ElementMapper& elementMapper_;
    const LookUpData lookUpData_;
    const LookUpCartesianData lookUpCartesianData_;
    const EclipseState& eclState_;
    std::vector<Scalar> thpresDefault_;
    std::vector<Scalar> thpres_;
    unsigned numEquilRegions_{};
    std::vector<unsigned short> elemEquilRegion_;

    // threshold pressure accross faults. EXPERIMENTAL!
    std::vector<Scalar> thpresftValues_;
    std::vector<int> cartElemFaultIdx_;

    bool enableThresholdPressure_;
};

} // namespace Opm

#endif
