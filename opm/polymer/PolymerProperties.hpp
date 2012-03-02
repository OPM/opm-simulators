/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_POLYMERPROPERTIES_HEADER_INCLUDED
#define OPM_POLYMERPROPERTIES_HEADER_INCLUDED


#include <vector>
#include <opm/core/utility/linearInterpolation.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>


namespace Opm
{

    class PolymerProperties
    {
    public:
        PolymerProperties()
        {
        }

        PolymerProperties(double c_max,
                          double mix_param,
                          double rock_density,
                          double dead_pore_vol,
                          const std::vector<double>& c_vals_visc,
                          const std::vector<double>& visc_mult_vals,
                          const std::vector<double>& c_vals_ads,
                          const std::vector<double>& ads_vals)
            : c_max_(c_max),
              mix_param_(mix_param),
              rock_density_(rock_density),
              dead_pore_vol_(dead_pore_vol),
              c_vals_visc_(c_vals_visc),
              visc_mult_vals_(visc_mult_vals),
              c_vals_ads_(c_vals_ads),
              ads_vals_(ads_vals)
        {
        }

        PolymerProperties(const EclipseGridParser& gridparser)
        {
            readFromDeck(gridparser);
        }

        void set(double c_max,
                 double mix_param,
                 double rock_density,
                 double dead_pore_vol,
                 const std::vector<double>& c_vals_visc,
                 const std::vector<double>& visc_mult_vals,
                 const std::vector<double>& c_vals_ads,
                 const std::vector<double>& ads_vals)
        {
            c_max_ = c_max;
            mix_param_ = mix_param;
            rock_density_ = rock_density;
            dead_pore_vol_ = dead_pore_vol;
            c_vals_visc_ = c_vals_visc;
            visc_mult_vals_ = visc_mult_vals;
            c_vals_ads_ = c_vals_ads;
            ads_vals_ = ads_vals;
        }

        void readFromDeck(const EclipseGridParser& gridparser)
        {

            // We assume NTMISC=1
            const std::vector<double>& plymax = gridparser.getPLYMAX().plymax_;
            c_max_ = plymax[0];
            const std::vector<double>& tlmixpar = gridparser.getTLMIXPAR().tlmixpar_;
            mix_param_ = tlmixpar[0];

            // We assume NTSFUN=1
            const std::vector<double>& plyrock = gridparser.getPLYROCK().plyrock_;
            dead_pore_vol_ = plyrock[0];
            rock_density_ = plyrock[3];

            // We assume NTPVT=1
            const PLYVISC& plyvisc = gridparser.getPLYVISC();
            c_vals_visc_ = plyvisc.concentration_;
            visc_mult_vals_ = plyvisc.factor_;

            // We assume NTSFUN=1
            const PLYADS& plyads = gridparser.getPLYADS();
            c_vals_ads_ = plyads.local_concentration_;
            ads_vals_ = plyads.adsorbed_concentration_;

        }

        double cMax() const
        {
            return c_max_;
        }

        double mixParam() const
        {
            return mix_param_;
        }

        double rockDensity() const
        {
            return rock_density_;
        };

        double deadPoreVol() const
        {
            return dead_pore_vol_;
        }

        double viscMult(double c) const
        {
            return Opm::linearInterpolation(c_vals_visc_, visc_mult_vals_, c);
        }

        double viscMultWithDer(double c, double* der) const
        {
            *der = Opm::linearInterpolationDerivative(c_vals_visc_, visc_mult_vals_, c);
            return Opm::linearInterpolation(c_vals_visc_, visc_mult_vals_, c);
        }

        double adsorbtion(double c) const
        {
            return Opm::linearInterpolation(c_vals_ads_, ads_vals_, c);
        }

        double adsorbtionWithDer(double c, double* der) const
        {
            *der = Opm::linearInterpolationDerivative(c_vals_ads_, ads_vals_, c);
            return Opm::linearInterpolation(c_vals_ads_, ads_vals_, c);
        }

    private:
        double c_max_;
        double mix_param_;
        double rock_density_;
        double dead_pore_vol_;
        std::vector<double> c_vals_visc_;
        std::vector<double> visc_mult_vals_;
        std::vector<double> c_vals_ads_;
        std::vector<double> ads_vals_;
    };

} // namespace Opm

#endif // OPM_POLYMERPROPERTIES_HEADER_INCLUDED
