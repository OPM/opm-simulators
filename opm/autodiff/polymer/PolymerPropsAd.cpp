
#include <cmath>
#include <vector>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/polymer/PolymerPropsAd.cpp>

namespace Opm {


    typedef PolymerPropsAd::ADB ADB;
    typedef PolymerPropsAd::V V;


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
    

    double PolymerProsAd::num_cells() const
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


    ADB PolymerPropsAd::effectiveInvWaterVisc(const ADB& c,
					      const double* visc)
    {
	const int n = c.size();
	V inv_mu_w_eff(n);
	V dinv_mu_w_eff(n);
	for (int i = 0; i < n; ++i) {
	    double im, dim;
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
    


} // namespace Opm
