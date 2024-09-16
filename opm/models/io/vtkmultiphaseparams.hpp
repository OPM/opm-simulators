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
 * \copydoc Opm::VtkMultiPhaseModule
 */
#ifndef OPM_VTK_MULTI_PHASE_PARAMS_HPP
#define OPM_VTK_MULTI_PHASE_PARAMS_HPP

namespace Opm::Parameters {

// set default values for what quantities to output
struct VtkWriteExtrusionFactor { static constexpr bool value = false; };
struct VtkWritePressures { static constexpr bool value = true; };
struct VtkWriteDensities { static constexpr bool value = true; };
struct VtkWriteSaturations { static constexpr bool value = true; };
struct VtkWriteMobilities { static constexpr bool value = false; };
struct VtkWriteRelativePermeabilities { static constexpr bool value = true; };
struct VtkWriteViscosities { static constexpr bool value = false; };
struct VtkWriteAverageMolarMasses { static constexpr bool value = false; };
struct VtkWritePorosity { static constexpr bool value = true; };
struct VtkWriteIntrinsicPermeabilities { static constexpr bool value = false; };
struct VtkWritePotentialGradients { static constexpr bool value = false; };
struct VtkWriteFilterVelocities { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

/*!
 * \brief Struct holding the parameters for VtkMultiPhaseModule.
 */
struct VtkMultiPhaseParams
{
    //! \brief Registers the parameters in parameter system.
    static void registerParameters();

    //! \brief Reads the parameter values from the parameter system.
    void read();

    bool extrusionFactorOutput_;
    bool pressureOutput_;
    bool densityOutput_;
    bool saturationOutput_;
    bool mobilityOutput_;
    bool relativePermeabilityOutput_;
    bool viscosityOutput_;
    bool averageMolarMassOutput_;
    bool porosityOutput_;
    bool intrinsicPermeabilityOutput_;
    bool velocityOutput_;
    bool potentialGradientOutput_;
};

} // namespace Opm

#endif // OPM_VTK_MULTI_PHASE_MODULE_HPP
