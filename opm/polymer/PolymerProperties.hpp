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


#include <cmath>
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
                          double res_factor,
                          double c_max_ads,
                          int ads_index,
                          const std::vector<double>& c_vals_visc,
                          const std::vector<double>& visc_mult_vals,
                          const std::vector<double>& c_vals_ads,
                          const std::vector<double>& ads_vals
                          )
            : c_max_(c_max),
              mix_param_(mix_param),
              rock_density_(rock_density),
              dead_pore_vol_(dead_pore_vol),
              res_factor_(res_factor),
              c_max_ads_(c_max_ads),
              ads_index_(ads_index),
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
                 double res_factor,
                 double c_max_ads,
                 int ads_index,
                 const std::vector<double>& c_vals_visc,
                 const std::vector<double>& visc_mult_vals,
                 const std::vector<double>& c_vals_ads,
                 const std::vector<double>& ads_vals
                 )
        {
            c_max_ = c_max;
            mix_param_ = mix_param;
            rock_density_ = rock_density;
            dead_pore_vol_ = dead_pore_vol;
            res_factor_ = res_factor;
            c_max_ads_ = c_max_ads;
            c_vals_visc_ = c_vals_visc;
            visc_mult_vals_ = visc_mult_vals;
            c_vals_ads_ = c_vals_ads;
            ads_vals_ = ads_vals;
            ads_index_ = ads_index;
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
            res_factor_ = plyrock[2];
            rock_density_ = plyrock[3];
            ads_index_ = plyrock[4];
            c_max_ads_ = plyrock[5];

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

        double resFactor() const
        {
            return res_factor_;
        }

        double cMaxAds() const
        {
            return c_max_ads_;
        }

        int adsIndex() const
        {
            return ads_index_;
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

        double simpleAdsorbtion(double c) const
        {
            return Opm::linearInterpolation(c_vals_ads_, ads_vals_, c);
        }

        double simpleAdsorbtionWithDer(double c, double* der) const
        {
            *der = Opm::linearInterpolationDerivative(c_vals_ads_, ads_vals_, c);
            return Opm::linearInterpolation(c_vals_ads_, ads_vals_, c);
        }

        double adsorbtion(double c, double cmax) const
        {
            if (ads_index_ == 1) {
                return simpleAdsorbtion(c);
            } else if (ads_index_ == 2) {
                return simpleAdsorbtion(std::max(c, cmax));
            } else {
                THROW("Invalid Adsoption index");
            }
        }

        double adsorbtionWithDer(double c, double cmax, double* der) const
        {
            if (ads_index_ == 1) {
                return simpleAdsorbtionWithDer(c, der);
            } else if (ads_index_ == 2) {
                if (c < cmax) {
                    *der = 0;
                    return simpleAdsorbtion(cmax);
                } else {
                    return simpleAdsorbtionWithDer(c, der);
                }
            } else {
                THROW("Invalid Adsoption index");
            }
        }


	void effectiveInvVisc(const double c, const double* visc, double* inv_visc_eff) const
	{
	    double cbar = c/c_max_;
	    double mu_w = visc[0];
	    double mu_m = viscMult(c)*mu_w;
	    double mu_p = viscMult(c_max_)*mu_w;
	    double mu_m_omega = std::pow(mu_m, mix_param_);
	    double mu_w_e   = mu_m_omega*std::pow(mu_w, 1.0 - mix_param_);
	    double mu_p_eff = mu_m_omega*std::pow(mu_p, 1.0 - mix_param_);
	    double inv_mu_w_eff = (1.0 - cbar)/mu_w_e + cbar/mu_p_eff;
	    inv_visc_eff[0] = inv_mu_w_eff;
	    inv_visc_eff[1] = 1.0/visc[1];
	}
        
        void effectiveMobilities(const double c,
                                 const double* visc,
                                 const double* relperm,
                                 const double* drelpermds,
                                 std::vector<double>& mob, 
                                 std::vector<double>& dmobds,
                                 double& dmobwatdc) const 
        {
            double cbar = c/c_max_;
            double mu_w = visc[0];
            double mu_m_dc; // derivative of mu_m with respect to c
            double mu_m = viscMultWithDer(c, &mu_m_dc)*mu_w;
            mu_m_dc *= mu_w;
            double mu_p = viscMult(c_max_)*mu_w;
            double omega = mix_param_;
            double mu_w_e   = std::pow(mu_m, omega)*std::pow(mu_w, 1 - omega);
            double mu_w_e_dc = omega*mu_m_dc*std::pow(mu_m, omega - 1)*std::pow(mu_w, 1 - omega);
            double mu_p_eff = std::pow(mu_m, omega)*std::pow(mu_p, 1 - omega);
            double mu_p_eff_dc = omega*mu_m_dc*std::pow(mu_m, omega - 1)*std::pow(mu_p, 1 - omega);
            double mu_w_eff = 1./((1 - cbar)/mu_w_e + cbar/mu_p_eff);
            double mu_w_eff_dc = -1./c_max_*mu_w_eff*mu_w_eff*(1./mu_p_eff - 1./mu_w_e)
                + (1-cbar)*(mu_w_eff*mu_w_eff/(mu_w_e*mu_w_e))*mu_w_e_dc
                + cbar*(mu_w_eff*mu_w_eff/(mu_p_eff*mu_p_eff))*mu_p_eff_dc;
            double visc_eff[2] = { mu_w_eff, visc[1] };

            dmobwatdc = - mob[0]*mu_w_eff_dc/(mu_w_eff*mu_w_eff);

            mob[0] = relperm[0]/visc_eff[0];
            mob[1] = relperm[1]/visc_eff[1];

            dmobds[0*2 + 0] = drelpermds[0*2 + 0]/visc_eff[0];
            dmobds[0*2 + 1] = drelpermds[0*2 + 1]/visc_eff[1];
            dmobds[1*2 + 0] = drelpermds[1*2 + 0]/visc_eff[0];
            dmobds[1*2 + 1] = drelpermds[1*2 + 1]/visc_eff[1];
        }

        void computeMc(const double& c,  
                       double& mc,
                       double& dmcdc) const 
        {
            double cbar = c/c_max_;
            double omega = mix_param_;
            double r = std::pow(viscMult(c_max_), 1 - omega); // viscMult(c_max_)=mu_p/mu_w
            mc = c/(cbar + (1 - cbar)*r);
            dmcdc = r/std::pow(cbar + (1 - cbar)*r, 2);
        }


    private:
        double c_max_;
        double mix_param_;
        double rock_density_;
        double dead_pore_vol_;
        double res_factor_;
        double c_max_ads_;
        int ads_index_;
        std::vector<double> c_vals_visc_;
        std::vector<double> visc_mult_vals_;
        std::vector<double> c_vals_ads_;
        std::vector<double> ads_vals_;
    };

} // namespace Opm

#endif // OPM_POLYMERPROPERTIES_HEADER_INCLUDED
