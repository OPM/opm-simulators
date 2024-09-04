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
 * \brief Contains the parameters required to extend the black-oil model by solvents.
 */
#ifndef OPM_BLACK_OIL_SOLVENT_PARAMS_HPP
#define OPM_BLACK_OIL_SOLVENT_PARAMS_HPP

#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/H2GasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BrineCo2Pvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BrineH2Pvt.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>

namespace Opm {

#if HAVE_ECL_INPUT
class EclipseState;
class Schedule;
#endif

//! \brief Struct holding the parameters for the BlackOilSolventModule class.
template<class Scalar>
struct BlackOilSolventParams
{
    using BrineCo2Pvt = ::Opm::BrineCo2Pvt<Scalar>;
    using BrineH2Pvt = ::Opm::BrineH2Pvt<Scalar>;
    using Co2GasPvt = ::Opm::Co2GasPvt<Scalar>;
    using H2GasPvt = ::Opm::H2GasPvt<Scalar>;
    using SolventPvt = ::Opm::SolventPvt<Scalar>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

#if HAVE_ECL_INPUT
    template<bool enableSolvent>
    void initFromState(const EclipseState& eclState, const Schedule& schedule);
#endif

    BrineCo2Pvt brineCo2Pvt_;
    BrineH2Pvt brineH2Pvt_;
    Co2GasPvt co2GasPvt_;
    H2GasPvt h2GasPvt_;
    SolventPvt solventPvt_;

    std::vector<TabulatedFunction> ssfnKrg_; // the krg(Fs) column of the SSFN table
    std::vector<TabulatedFunction> ssfnKrs_; // the krs(Fs) column of the SSFN table
    std::vector<TabulatedFunction> sof2Krn_; // the krn(Sn) column of the SOF2 table
    std::vector<TabulatedFunction> misc_;    // the misc(Ss) column of the MISC table
    std::vector<TabulatedFunction> pmisc_;   // the pmisc(pg) column of the PMISC table
    std::vector<TabulatedFunction> msfnKrsg_; // the krsg(Ssg) column of the MSFN table
    std::vector<TabulatedFunction> msfnKro_; // the kro(Ssg) column of the MSFN table
    std::vector<TabulatedFunction> sorwmis_; // the sorwmis(Sw) column of the SORWMIS table
    std::vector<TabulatedFunction> sgcwmis_; // the sgcwmis(Sw) column of the SGCWMIS table

    std::vector<Scalar> tlMixParamViscosity_; // Todd-Longstaff mixing parameter for viscosity
    std::vector<Scalar> tlMixParamDensity_;   //  Todd-Longstaff mixing parameter for density
    std::vector<TabulatedFunction> tlPMixTable_; // the tlpmixpa(Po) column of the TLPMIXPA table

    bool isMiscible_;
    bool rsSolw_active_ = false;
    bool co2sol_;
    bool h2sol_;

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the SSFN of any region.
     */
    void setNumSatRegions(unsigned numRegions);

    /*!
     * \brief Specify miscible relative permeability multipliers of a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    void setMsfn(unsigned satRegionIdx,
                 const TabulatedFunction& msfnKrsg,
                 const TabulatedFunction& msfnKro);
};

} // namespace Opm

#endif // OPM_BLACK_OIL_SOLVENT_PARAMS_HPP
