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

#include <opm/core/props/IncompPropertiesSinglePhase.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/Units.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/parser/eclipse/Deck/DeckRecord.hpp>

namespace Opm
{
    IncompPropertiesSinglePhase::IncompPropertiesSinglePhase(Opm::DeckConstPtr deck,
                                                             Opm::EclipseStateConstPtr eclState,
                                                             const UnstructuredGrid& grid)
    {
        rock_.init(eclState, grid.number_of_cells, grid.global_cell, grid.cartdims);

        if (deck->hasKeyword("DENSITY")) {
            Opm::DeckRecordConstPtr densityRecord = deck->getKeyword("DENSITY")->getRecord(0);
            surface_density_ = densityRecord->getItem("OIL")->getSIDouble(0);
        } else {
            surface_density_ = 1000.0;
            OPM_MESSAGE("Input is missing DENSITY -- using a standard density of "
                        << surface_density_ << ".\n");
        }

        // This will be modified if we have a PVCDO specification.
        reservoir_density_ = surface_density_;

        if (deck->hasKeyword("PVCDO")) {
            Opm::DeckRecordConstPtr pvcdoRecord = deck->getKeyword("PVCDO")->getRecord(0);
            if (pvcdoRecord->getItem("OIL_COMPRESSIBILITY")->getSIDouble(0) != 0.0 ||
                pvcdoRecord->getItem("OIL_VISCOSIBILITY")->getSIDouble(0) != 0.0) {
                OPM_MESSAGE("Compressibility effects in PVCDO are ignored.");
            }
            reservoir_density_ /= pvcdoRecord->getItem("OIL_VOL_FACTOR")->getSIDouble(0);
            viscosity_ = pvcdoRecord->getItem("OIL_VISCOSITY")->getSIDouble(0);
        } else {
            viscosity_ = 1.0 * prefix::centi*unit::Poise;
            OPM_MESSAGE("Input is missing PVCDO -- using a standard viscosity of "
                        << viscosity_ << " and reservoir density equal to surface density.\n");
        }
    }

    IncompPropertiesSinglePhase::~IncompPropertiesSinglePhase()
    {
    }


    /// \return   D, the number of spatial dimensions.
    int IncompPropertiesSinglePhase::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int IncompPropertiesSinglePhase::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* IncompPropertiesSinglePhase::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* IncompPropertiesSinglePhase::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int IncompPropertiesSinglePhase::numPhases() const
    {
        return 1;
    }

    /// \return Array of P viscosity values.
    const double* IncompPropertiesSinglePhase::viscosity() const
    {
        return &viscosity_;
    }

    /// \return Array of P density values.
    const double* IncompPropertiesSinglePhase::density() const
    {
        return &reservoir_density_;
    }

    /// \return Array of P density values.
    const double* IncompPropertiesSinglePhase::surfaceDensity() const
    {
        return &surface_density_;
    }

    /// Relative permeability. Always returns 1 (and 0 for derivatives).
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of n saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of n relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of n relperm derivative values,
    ///                    array must be valid before calling.
    void IncompPropertiesSinglePhase::relperm(const int n,
                                              const double* /* s */,
                                              const int* /* cells */,
                                              double* kr,
                                              double* dkrds) const
    {
        std::fill(kr, kr + n, 1.0);
        if (dkrds) {
            std::fill(dkrds, dkrds + n, 0.0);
        }
    }


    /// Capillary pressure. Always returns zero.
    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of n saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of n capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of n derivative values,
    ///                    array must be valid before calling.
    void IncompPropertiesSinglePhase::capPress(const int n,
                                               const double* /* s */,
                                               const int* /* cells */,
                                               double* pc,
                                               double* dpcds) const
    {
        std::fill(pc, pc + n, 0.0);
        if (dpcds) {
            std::fill(dpcds, dpcds + n, 0.0);
        }
    }


    /// Obtain the range of allowable saturation values.
    /// Saturation range is just the point 1 for this class
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of n minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of n maximum s values, array must be valid before calling.
    void IncompPropertiesSinglePhase::satRange(const int n,
                                               const int* /* cells */,
                                               double* smin,
                                               double* smax) const
    {
        std::fill(smin, smin + n, 1.0);
        std::fill(smax, smax + n, 1.0);
    }

} // namespace Opm

