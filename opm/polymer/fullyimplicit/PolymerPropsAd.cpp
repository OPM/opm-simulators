/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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
#include <cmath>
#include <vector>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>

namespace Opm {

    typedef PolymerPropsAd::ADB ADB;
    typedef PolymerPropsAd::V V;





    double
    PolymerPropsAd::rockDensity() const
    {
        return polymer_props_.rockDensity();
    }





    double
    PolymerPropsAd::deadPoreVol() const
    {
        return polymer_props_.deadPoreVol();
    }





	double
	PolymerPropsAd::cMax() const
	{
		return polymer_props_.cMax();
	}

    const std::vector<double>&
    PolymerPropsAd::shearWaterVelocity() const
    {
        return polymer_props_.shearWaterVelocity();
    }

    const std::vector<double>&
    PolymerPropsAd::shearViscosityReductionFactor() const
    {
        return polymer_props_.shearViscosityReductionFactor();
    }

    double
    PolymerPropsAd::plyshlogRefConc() const
    {
        return polymer_props_.plyshlogRefConc();
    }

    bool
    PolymerPropsAd::hasPlyshlogRefSalinity() const
    {
        return polymer_props_.hasPlyshlogRefSalinity();
    }

    bool
    PolymerPropsAd::hasPlyshlogRefTemp() const
    {
        return polymer_props_.hasPlyshlogRefTemp();
    }

    double
    PolymerPropsAd::plyshlogRefSalinity() const
    {
        return polymer_props_.plyshlogRefSalinity();
    }

    double
    PolymerPropsAd::plyshlogRefTemp() const
    {
        return polymer_props_.plyshlogRefTemp();
    }

    double
    PolymerPropsAd::viscMult(double c) const
    {
        return polymer_props_.viscMult(c);
    }

    V
    PolymerPropsAd::viscMult(const V& c) const
    {
        int nc = c.size();
        V visc_mult(nc);
        for (int i = 0; i < nc; ++i) {
            visc_mult[i] = polymer_props_.viscMult(c[i]);
        }
        return visc_mult;
    }


	

    PolymerPropsAd::PolymerPropsAd(const PolymerProperties& polymer_props)
        : polymer_props_ (polymer_props)
    {
    }





    PolymerPropsAd::~PolymerPropsAd()
    {
    }





    V PolymerPropsAd::effectiveInvWaterVisc(const V& c,
                                            const double* visc) const
    {
        const int nc = c.size();
        V inv_mu_w_eff(nc);
        for (int i = 0; i < nc; ++i) {
            double im = 0;
            polymer_props_.effectiveInvVisc(c(i), visc, im);
            inv_mu_w_eff(i) = im;
        }

        return inv_mu_w_eff;
    }





    ADB PolymerPropsAd::effectiveInvWaterVisc(const ADB& c,
	                    				      const double* visc) const
    {
	    const int nc = c.size();
    	V inv_mu_w_eff(nc);
    	V dinv_mu_w_eff(nc);
    	for (int i = 0; i < nc; ++i) {
    	    double im = 0, dim = 0;
    	    polymer_props_.effectiveInvViscWithDer(c.value()(i), visc, im, dim);
    	    inv_mu_w_eff(i) = im;
    	    dinv_mu_w_eff(i) = dim;
    	}
        ADB::M dim_diag = spdiag(dinv_mu_w_eff);
        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dim_diag * c.derivative()[block];
        }
        return ADB::function(std::move(inv_mu_w_eff), std::move(jacs));
    }





    V PolymerPropsAd::polymerWaterVelocityRatio(const V& c) const
    {
        const int nc = c.size();
        V mc(nc);

        for (int i = 0; i < nc; ++i) {
            double m = 0;
            polymer_props_.computeMc(c(i), m);
            mc(i) = m;
        }
       
       return mc;
    }





    ADB PolymerPropsAd::polymerWaterVelocityRatio(const ADB& c) const
    {
    
        const int nc = c.size();
        V mc(nc);
        V dmc(nc);
        
        for (int i = 0; i < nc; ++i) {
            double m = 0;
            double dm = 0;
            polymer_props_.computeMcWithDer(c.value()(i), m, dm);

            mc(i) = m;
            dmc(i) = dm;
        }

        ADB::M dmc_diag = spdiag(dmc);
        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmc_diag * c.derivative()[block];
        }

        return ADB::function(std::move(mc), std::move(jacs));
    }





    V PolymerPropsAd::adsorption(const V& c, const V& cmax_cells) const
    {
        const int nc = c.size();
        V ads(nc);

        for (int i = 0; i < nc; ++i) {
            double c_ads = 0;
            polymer_props_.adsorption(c(i), cmax_cells(i), c_ads);

            ads(i) = c_ads;
        }

        return ads;
    }





    ADB PolymerPropsAd::adsorption(const ADB& c, const ADB& cmax_cells) const
    {
        const int nc = c.value().size();

        V ads(nc);
        V dads(nc);

        for (int i = 0; i < nc; ++i) {
            double c_ads = 0;
            double dc_ads = 0;
            polymer_props_.adsorptionWithDer(c.value()(i), cmax_cells.value()(i), c_ads, dc_ads);
            ads(i) = c_ads;
            dads(i) = dc_ads;
        }

        ADB::M dads_diag = spdiag(dads);
        int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dads_diag * c.derivative()[block];
        }

        return ADB::function(std::move(ads), std::move(jacs));
    }





    V
    PolymerPropsAd::effectiveRelPerm(const V& c, 
                                     const V& cmax_cells,
                                     const V& krw) const
    {
        const int nc = c.size();

        V one  = V::Ones(nc);
        V ads = adsorption(c, cmax_cells);
        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor -1.) / max_ads;
        V rk = one + factor * ads;

        return krw / rk;
    }





    ADB
    PolymerPropsAd::effectiveRelPerm(const ADB& c,
                                     const ADB& cmax_cells,
                                     const ADB& krw) const
    {
        const int nc = c.value().size();
        V one = V::Ones(nc);
        ADB ads = adsorption(c, cmax_cells);
        V krw_eff = effectiveRelPerm(c.value(), cmax_cells.value(), krw.value());

        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor - 1.) / max_ads;
        ADB rk = one + ads * factor; 
        
        return krw / rk;
    }


    bool
    PolymerPropsAd::computeShearMultLog(std::vector<double>& water_vel, std::vector<double>& visc_mult, std::vector<double>& shear_mult) const
    {
        return polymer_props_.computeShearMultLog(water_vel, visc_mult, shear_mult);
    }


}// namespace Opm
