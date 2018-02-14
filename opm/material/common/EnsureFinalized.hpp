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
 * \copydoc Opm::EnsureFinalized
 */
#ifndef OPM_MATERIAL_ENSURE_FINALIZED_HPP
#define OPM_MATERIAL_ENSURE_FINALIZED_HPP

#include <cassert>
#include <stdexcept>

// TODO: move this variable to config.h
#define OPM_CHECK_PARAM_FINALIZED 1

#if ! defined(NDEBUG) && OPM_CHECK_PARAM_FINALIZED
#define USE_OPM_CHECK_PARAM_FINALIZED 1
#endif

namespace Opm {

/*!
 * \brief Default implementation for asserting finalization of parameter objects.
 *
 */
class EnsureFinalized
{
#if USE_OPM_CHECK_PARAM_FINALIZED
    bool finalized_;
#endif

protected:
    /*!
     * \brief The default constructor.
     */
    EnsureFinalized()
#if USE_OPM_CHECK_PARAM_FINALIZED
        : finalized_( false )
#endif
    {
    }

    void check() const
    {
#if USE_OPM_CHECK_PARAM_FINALIZED
        if (!finalized_)
            throw std::runtime_error("Parameter class has not been finalized before usage!");
#endif
    }

public:
    /*!
     * \brief Mark the object as finalized.
     */
    void finalize()
    {
#if USE_OPM_CHECK_PARAM_FINALIZED
        finalized_ = true;
#endif
    }
};

#undef USE_OPM_CHECK_PARAM_FINALIZED

} // namespace Opm
#endif
