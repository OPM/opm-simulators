
#include <cmath>
#include <vector>
#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/polymer/fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>

namespace Opm {


    typedef PolymerPropsAd::ADB ADB;
    typedef PolymerPropsAd::V V;

/*
    PolymerPropsAd::PolymerPropsAd()
    {
    }


    PolymerPropsAd::PolymerPropsAd(const int num_cells,
                                   const double c_max,
                                   const double mix_param,
                                   const double std::vector<double>& c_vals_visc,
                                   const double std::vector<double>& visc_mult_vals)
        : num_cells_ (num_cells)
        , c_max_ (c_max)
        , mix_param_(mix_param)
        , c_vals_visc_ (c_vals_visc)
        , visc_mult_vals_ (visc_mult_vals)
    {
    }
    

    double PolymerPropsAd::num_cells() const
    {
        return num__cells_;
    }




    double PolymerPropsAd::cMax() const
    {
        return c_max_;
    }



    double PolymerPropsAd::mixParam() const
    {
        return mix_param_;
    }



    V PolymerPropsAd::muM(const V& c,
                          const double* visc) const
    {
        const int nc = num_cells();
        assert(nc == c.size());
        std::vector<double> mu(nc);

        for (int i = 0; i < nc; ++i) {
            mu[i] = Opm::linearInterpolation(c_vals_visc_, visc_mult_vals_, c(i));
        }

        const V muM = Eigen::Map<const V>(&mu[0], nc);
        
        const double mu_w = visc[0];

        return muM * mu_w;
    }

   

    ADB PolymerPropsAd::muM(const ADB& c,
                            const double* visc) const
    {
        const int nc = num_cells();
        assert(nc == c.size());
        
        V mu_m = muM(c.value());

        std::vector<double> dmu_dc(nc);

        for (int i = 0; i < nc; ++i) {
            dmu_dc[i] = Opm::linearInterpolationDerivative(c_vals_visc_, visc_mult_vals_, c.value()(i));
        }      
        
        const V dmu = Eigen::Map<const V>(&dmu_dc[0], nc);

        ADB::M dmu_diag = spdiag(dmu);
        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);

        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmu_diag * c.derivative()[block];
        }
        
        const double mu_w = visc[0]

        return ADB::function(mu_m, jacs) * mu_w;
    }



    V PolymerPropsAd::ToddLongstaff(const double mix_param,
                                    const V& muM,
                                    const V& mu) const
    {
        const int nc = num_cells();
        assert(nc == muM.size());

        std::vector<double> mueff(nc);
        const double omega = mix_param;

        for (int i = 0; i < nc; ++i) {
            mueff[i] = std::pow(muM(i),omega) * std::pow(mu(i), 1. - omega);
        }

        const V muEff = Eigen::Map<const V>(&mueff[0], nc);

        return muEff;
    }


    ADB PolymerPropsAd::ToddLongstaff(const double mix_param,
                                      const ADB& muM,
                                      const ADB& mu) const
    {
        const int nc = num_cells();

    }


    V PolymerPropsAd::muPolyEff(const double mix_param,
                                const V& muM,
                                const V& muPoly) const
    {
    
        return ToddLongstaff(mix_param, muM, muPoly);
    }

    V PolymerPropsAd::muWatEff(const double mix_param,
                               const std::vecotr<double>& c_max,
                               const V& c,
                               const V& muM,
                               const V& muWat,
                               const V& muPolyEff) const
    {
        const int nc = num_cells();
        assert(nc == c.size());
        V muWate = ToddLongstaff(mix_param, muM, muWat);

//        V cmax = V::Constant(nc, 1, c_max);
        const V cmax = Eigen::Map<const V>(&c_max[0], nc);
        const V one = V::Ones(nc, 1);
        V inv_muWatEff = (one - c / cmax) / muWate + c / cmax / muPolyEff;

        V muWatEff = one / inv_muWatEff;

        return muWatEff;

    }
*/
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
        return ADB::function(inv_mu_w_eff, jacs);
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

        return ADB::function(mc, jacs);
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

        return ADB::function(ads, jacs);
    }


    V
    PolymerPropsAd::effectiveRelPerm(const V& c, 
                                     const V& cmax_cells,
                                     const V& krw) const
    {
        const int nc = c.size();

        V one  = V::Ones(nc);
        V ads = adsorption(c);
        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor -1.) / max_ads;
        V rk = one + factor * ads;
        V krw_eff = krw / rk;

        return eff_relperm;
    }


    ADB
    PolymerPropsAd::effectiveRelPerm(const ADB& c,
                                     const ADB& cmax_cells,
                                     const ADB& krw,
                                     const ADB& sw) const
    {
        const int nc = c.value().size();

        V one = V::Ones(nc);

        ADB ads = adsorption(c);
        V krw_eff = effectiveRelPerm(c.value(), cmax_cells.value(), krw.value());

        double max_ads = polymer_props_.cMaxAds();
        double res_factor = polymer_props_.resFactor();
        double factor = (res_factor - 1.) / max_ads;
        ADB rk = one + ads * factor; 
        ADB::M dkrw_ds = krw.derivative() / rk.derivative();
        ADB::M dkrw_dc = -krw.value() / (rk.value * rk.value()) * ads.derivative() * factor;

        const int num_blocks = c.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jac[block] = dkrw_ds * sw.derivative()[block] + dkrw_dc * c.derivative()[block];
        }

        return ADB::function(krw_eff, jacs);
    }

}// namespace Opm
