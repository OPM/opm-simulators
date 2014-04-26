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

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Utility/PlymaxTable.hpp>
#include <opm/parser/eclipse/Utility/TlmixparTable.hpp>
#include <opm/parser/eclipse/Utility/PlyrockTable.hpp>
#include <opm/parser/eclipse/Utility/PlyviscTable.hpp>
#include <opm/parser/eclipse/Utility/PlyadsTable.hpp>

#include <cmath>
#include <vector>


namespace Opm
{

    class PolymerProperties
    {
    public:
        PolymerProperties()
        {
        }

        enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

        /// Construct from parameters
        /// \param[in] c_max          Maximum polymer concentration used in computation of effective viscosity
        /// \param[in] mix_param      Mixing parameter
        /// \param[in] rock_density   Rock density
        /// \param[in] dead_pore_vol  Dead pore volume
        /// \param[in] res_factor     Residual resistance factor
        /// \param[in] c_max_ads      Maximum polymer adsorption value  used in computation of  the resistance factor
        /// \param[in] c_vals_visc    Array of concentration for effective vicosity multiplier
        /// \param[in] visc_mult_vals Array of effective vicosity multiplier
        /// \param[in] c_vals_ads     Array of concentration for adsorption values
        /// \param[in] ads_vals       Array of adsorption values
        /// \param[in] water_vel_vals_ Array of water phase velocity for shear
        /// \param[in] shear_vrf_vals_ Array of viscosity reduction factor
        PolymerProperties(double c_max,
                          double mix_param,
                          double rock_density,
                          double dead_pore_vol,
                          double res_factor,
                          double c_max_ads,
                          AdsorptionBehaviour ads_index,
                          const std::vector<double>& c_vals_visc,
                          const std::vector<double>& visc_mult_vals,
                          const std::vector<double>& c_vals_ads,
                          const std::vector<double>& ads_vals,
                          const std::vector<double>& water_vel_vals,
                          const std::vector<double>& shear_vrf_vals
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
              ads_vals_(ads_vals),
              water_vel_vals_(water_vel_vals),
              shear_vrf_vals_(shear_vrf_vals)
        {
        }

        PolymerProperties(Opm::DeckConstPtr deck)
        {
            readFromDeck(deck);
        }

        void set(double c_max,
                 double mix_param,
                 double rock_density,
                 double dead_pore_vol,
                 double res_factor,
                 double c_max_ads,
                 AdsorptionBehaviour ads_index,
                 const std::vector<double>& c_vals_visc,
                 const std::vector<double>& visc_mult_vals,
                 const std::vector<double>& c_vals_ads,
                 const std::vector<double>& ads_vals,
                 const std::vector<double>& water_vel_vals,
                 const std::vector<double>& shear_vrf_vals
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
            water_vel_vals_ = water_vel_vals;
            shear_vrf_vals_ = shear_vrf_vals;
        }

        void readFromDeck(Opm::DeckConstPtr deck)
        {
            // We assume NTMISC=1
            Opm::PlymaxTable plymaxTable(deck->getKeyword("PLYMAX"), /*tableIdx=*/0);
            Opm::TlmixparTable tlmixparTable(deck->getKeyword("TLMIXPAR"), /*tableIdx=*/0);

            // We also assume that each table has exactly one row...
            assert(plymaxTable.numRows() == 1);
            assert(tlmixparTable.numRows() == 1);

            c_max_ = plymaxTable.getPolymerConcentrationColumn()[0];
            mix_param_ = tlmixparTable.getViscosityParameterColumn()[0];

            // We assume NTSFUN=1
            Opm::PlyrockTable plyrockTable(deck->getKeyword("PLYROCK"), /*tableIdx=*/0);

            // We also assume that each table has exactly one row...
            assert(plyrockTable.numRows() == 1);

            dead_pore_vol_ = plyrockTable.getDeadPoreVolumeColumn()[0];
            res_factor_ = plyrockTable.getResidualResistanceFactorColumn()[0];
            rock_density_ = plyrockTable.getRockDensityFactorColumn()[0];
            ads_index_ = static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[0]);
            c_max_ads_ = plyrockTable.getMaxAdsorbtionColumn()[0];

            // We assume NTPVT=1
            Opm::PlyviscTable plyviscTable(deck->getKeyword("PLYVISC"), /*tableIdx=*/0);

            // We also assume that each table has exactly one row...
            assert(plyviscTable.numRows() == 1);

            c_vals_visc_[0] = plyviscTable.getPolymerConcentrationColumn()[0];
            visc_mult_vals_[0] =  plyviscTable.getViscosityMultiplierColumn()[0];

            // We assume NTSFUN=1
            Opm::PlyadsTable plyadsTable(deck->getKeyword("PLYADS"), /*tableIdx=*/0);

            // We also assume that each table has exactly one row...
            assert(plyadsTable.numRows() == 1);

            c_vals_ads_[0] = plyadsTable.getPolymerConcentrationColumn()[0];
            ads_vals_[0] = plyadsTable.getAdsorbedPolymerColumn()[0];
        }

        double cMax() const;

        double mixParam() const;

        double rockDensity() const;

        double deadPoreVol() const;

        double resFactor() const;

        double cMaxAds() const;

        int adsIndex() const;
        
        const std::vector<double>& shearWaterVelocity() const;

        double shearVrf(const double velocity) const;

        double shearVrfWithDer(const double velocity, double& der) const;

        double viscMult(double c) const;

        double viscMultWithDer(double c, double* der) const;

        void simpleAdsorption(double c, double& c_ads) const;

        void simpleAdsorptionWithDer(double c, double& c_ads,
                                     double& dc_ads_dc) const;

        void adsorption(double c, double cmax, double& c_ads) const;

        void adsorptionWithDer(double c, double cmax,
                               double& c_ads, double& dc_ads_dc) const;

        void effectiveVisc(const double c, const double* visc,
                                              double& mu_w_eff) const;

        void effectiveViscWithDer(const double c, const double* visc
                                                     , double& mu_w_eff
                                                     , double dmu_w_eff_dc) const;

        void effectiveInvVisc(const double c, const double* visc,
                                                 double& inv_mu_w_eff) const;

        void effectiveInvViscWithDer(const double c,
                                                        const double* visc,
                                                        double& inv_mu_w_eff,
                                                        double& dinv_mu_w_eff_dc) const;
        void effectiveRelperm(const double c,
                              const double cmax,
                              const double* relperm,
                              double& eff_relperm_wat) const;

        void effectiveRelpermWithDer (const double c,
                                      const double cmax,
                                      const double* relperm,
                                      const double* drelperm_ds,
                                      double& eff_relperm_wat,
                                      double& deff_relperm_wat_ds,
                                      double& deff_relperm_wat_dc) const;

        void effectiveMobilities(const double c,
                                 const double cmax,
                                 const double* visc,
                                 const double* relperm,
                                 double* mob) const;

        void effectiveMobilitiesWithDer(const double c,
                                        const double cmax,
                                        const double* visc,
                                        const double* relperm,
                                        const double* drelpermds,
                                        double* mob,
                                        double*  dmob_ds,
                                        double& dmobwatdc) const;

        void effectiveMobilitiesBoth(const double c,
                                     const double cmax,
                                     const double* visc,
                                     const double* relperm,
                                     const double* drelperm_ds,
                                     double* mob,
                                     double* dmob_ds,
                                     double& dmobwat_dc,
                                     bool if_with_der) const;

        void effectiveTotalMobility(const double c,
                                    const double cmax,
                                    const double* visc,
                                    const double* relperm,
                                    double& totmob) const;

        void effectiveTotalMobilityWithDer(const double c,
                                           const double cmax,
                                           const double* visc,
                                           const double* relperm,
                                           const double* drelpermds,
                                           double& totmob,
                                           double* dtotmob_dsdc) const;

        void effectiveTotalMobilityBoth(const double c,
                                        const double cmax,
                                        const double* visc,
                                        const double* relperm,
                                        const double* drelperm_ds,
                                        double& totmob,
                                        double* dtotmob_dsdc,
                                        bool if_with_der) const;

        void computeMc(const double& c, double& mc) const;

        void computeMcWithDer(const double& c, double& mc,
                              double& dmc_dc) const;

        void computeMcBoth(const double& c, double& mc,
                           double& dmc_dc, bool if_with_der) const;

    private:
        double c_max_;
        double mix_param_;
        double rock_density_;
        double dead_pore_vol_;
        double res_factor_;
        double c_max_ads_;
        AdsorptionBehaviour ads_index_;
        std::vector<double> c_vals_visc_;
        std::vector<double> visc_mult_vals_;
        std::vector<double> c_vals_ads_;
        std::vector<double> ads_vals_;
        std::vector<double> water_vel_vals_;
        std::vector<double> shear_vrf_vals_;
        void simpleAdsorptionBoth(double c, double& c_ads,
                                  double& dc_ads_dc, bool if_with_der) const;
        void adsorptionBoth(double c, double cmax,
                            double& c_ads, double& dc_ads_dc,
                            bool if_with_der) const;
        void effectiveInvViscBoth(const double c, const double* visc,
                                  double& inv_mu_w_eff,
                                  double& dinv_mu_w_eff_dc, bool if_with_der) const;
        void effectiveRelpermBoth(const double c,
                                  const double cmax,
                                  const double* relperm,
                                  const double* drelperm_ds,
                                  double& eff_relperm_wat,
                                  double& deff_relperm_wat_ds,
                                  double& deff_relperm_wat_dc,
                                  bool if_with_der) const;
    };

} // namespace Opm

#endif // OPM_POLYMERPROPERTIES_HEADER_INCLUDED
