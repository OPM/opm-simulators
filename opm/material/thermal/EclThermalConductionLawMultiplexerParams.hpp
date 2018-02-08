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
 * \copydoc Opm::EclThermalConductionLawMultiplexerParams
 */
#ifndef OPM_ECL_THERMAL_CONDUCTION_LAW_MULTIPLEXER_PARAMS_HPP
#define OPM_ECL_THERMAL_CONDUCTION_LAW_MULTIPLEXER_PARAMS_HPP

#include "EclThconrLawParams.hpp"
#include "EclThcLawParams.hpp"

#include <opm/material/common/EnsureFinalized.hpp>

#include <memory>

namespace Opm {

/*!
 * \brief The default implementation of a parameter object for the
 *        ECL thermal law.
 */
template <class ScalarT>
class EclThermalConductionLawMultiplexerParams : public EnsureFinalized
{
    typedef void* ParamPointerType;

public:
    typedef ScalarT Scalar;

    enum ThermalConductionApproach {
        undefinedApproach,
        thconrApproach, // keywords: THCONR, THCONSF
        thcApproach, // keywords: THCROCK, THCOIL, THCGAS, THCWATER
        nullApproach, // (no keywords)
    };

    typedef Opm::EclThconrLawParams<ScalarT> ThconrLawParams;
    typedef Opm::EclThcLawParams<ScalarT> ThcLawParams;

    EclThermalConductionLawMultiplexerParams(const EclThermalConductionLawMultiplexerParams&) = default;

    EclThermalConductionLawMultiplexerParams()
    { thermalConductionApproach_ = undefinedApproach; }

    ~EclThermalConductionLawMultiplexerParams()
    { destroy_(); }

    void setThermalConductionApproach(ThermalConductionApproach newApproach)
    {
        destroy_();

        thermalConductionApproach_ = newApproach;
        switch (thermalConductionApproach()) {
        case undefinedApproach:
            throw std::logic_error("Cannot set the approach for thermal conduction to 'undefined'!");

        case thconrApproach:
            realParams_ = new ThconrLawParams;
            break;

        case thcApproach:
            realParams_ = new ThcLawParams;
            break;

        case nullApproach:
            realParams_ = nullptr;
            break;
        }
    }

    ThermalConductionApproach thermalConductionApproach() const
    { return thermalConductionApproach_; }

    // get the parameter object for the THCONR case
    template <ThermalConductionApproach approachV>
    typename std::enable_if<approachV == thconrApproach, ThconrLawParams>::type&
    getRealParams()
    {
        assert(thermalConductionApproach() == approachV);
        return *static_cast<ThconrLawParams*>(realParams_);
    }

    template <ThermalConductionApproach approachV>
    typename std::enable_if<approachV == thconrApproach, const ThconrLawParams>::type&
    getRealParams() const
    {
        assert(thermalConductionApproach() == approachV);
        return *static_cast<const ThconrLawParams*>(realParams_);
    }

    // get the parameter object for the THC* case
    template <ThermalConductionApproach approachV>
    typename std::enable_if<approachV == thcApproach, ThcLawParams>::type&
    getRealParams()
    {
        assert(thermalConductionApproach() == approachV);
        return *static_cast<ThcLawParams*>(realParams_);
    }

    template <ThermalConductionApproach approachV>
    typename std::enable_if<approachV == thcApproach, const ThcLawParams>::type&
    getRealParams() const
    {
        assert(thermalConductionApproach() == approachV);
        return *static_cast<const ThcLawParams*>(realParams_);
    }

private:
    void destroy_()
    {
        switch (thermalConductionApproach()) {
        case undefinedApproach:
            break;

        case thconrApproach:
            delete static_cast<ThconrLawParams*>(realParams_);
            break;

        case thcApproach:
            delete static_cast<ThcLawParams*>(realParams_);
            break;

        case nullApproach:
            break;
        }

        thermalConductionApproach_ = undefinedApproach;
    }

    ThermalConductionApproach thermalConductionApproach_;
    ParamPointerType realParams_;
};

} // namespace Opm

#endif
