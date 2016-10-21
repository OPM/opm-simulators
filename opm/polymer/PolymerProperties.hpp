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
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyadsTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymaxTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyrockTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyviscTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/Units/Dimension.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>


#include <cmath>
#include <vector>
#include <opm/common/ErrorMacros.hpp>


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

        PolymerProperties(const Opm::Deck& deck, const Opm::EclipseState& eclipseState)
        {
            readFromDeck(deck, eclipseState);
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

        void readFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclipseState)
        {
            // We assume NTMISC=1
            const auto& tables = eclipseState.getTableManager();
            const auto& plymaxTable = tables.getPlymaxTables().getTable<PlymaxTable>(0);
            const auto& plmixparRecord = deck.getKeyword("PLMIXPAR").getRecord(0);

            // We also assume that each table has exactly one row...
            assert(plymaxTable.numRows() == 1);

            c_max_ = plymaxTable.getPolymerConcentrationColumn()[0];
            mix_param_ = plmixparRecord.getItem("TODD_LONGSTAFF").getSIDouble(0);

            // We assume NTSFUN=1
            const auto& plyrockTable = tables.getPlyrockTables().getTable<PlyrockTable>(0);

            // We also assume that each table has exactly one row...
            assert(plyrockTable.numRows() == 1);

            dead_pore_vol_ = plyrockTable.getDeadPoreVolumeColumn()[0];
            res_factor_ = plyrockTable.getResidualResistanceFactorColumn()[0];
            rock_density_ = plyrockTable.getRockDensityFactorColumn()[0];
            ads_index_ = static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[0]);
            c_max_ads_ = plyrockTable.getMaxAdsorbtionColumn()[0];

            // We assume NTPVT=1
            const auto& plyviscTable = tables.getPlyviscTables().getTable<PlyviscTable>(0);


            c_vals_visc_ = plyviscTable.getPolymerConcentrationColumn().vectorCopy( );
            visc_mult_vals_ =  plyviscTable.getViscosityMultiplierColumn().vectorCopy( );

            // We assume NTSFUN=1
            const auto& plyadsTable = tables.getPlyadsTables().getTable<PlyadsTable>(0);

            c_vals_ads_ = plyadsTable.getPolymerConcentrationColumn().vectorCopy( );
            ads_vals_ = plyadsTable.getAdsorbedPolymerColumn().vectorCopy( );

            has_plyshlog_ = deck.hasKeyword("PLYSHLOG");
            has_shrate_ = deck.hasKeyword("SHRATE");

            if (has_plyshlog_) {
                // Assuming NTPVT == 1 always
                const auto& plyshlogTable = tables.getPlyshlogTables().getTable<PlyshlogTable>(0);

                water_vel_vals_ = plyshlogTable.getWaterVelocityColumn().vectorCopy( );
                shear_vrf_vals_ = plyshlogTable.getShearMultiplierColumn().vectorCopy( );

                // do the unit version here for the water_vel_vals_
                Opm::UnitSystem unitSystem = deck.getActiveUnitSystem();
                double siFactor;
                if (has_shrate_) {
                    siFactor = unitSystem.parse("1/Time").getSIScaling();
                    const auto& shrateKeyword = deck.getKeyword("SHRATE");
                    std::vector<double> shrate_readin = shrateKeyword.getSIDoubleData();
                    if (shrate_readin.size() == 1) {
                        shrate_ = shrate_readin[0];
                    } else if (shrate_readin.size() == 0) {
                        shrate_ = 4.8; // default value
                    } else {
                        OPM_THROW(std::logic_error, "Only NTPVT == 1 is allowed for SHRATE keyword now !\n");
                    }
                } else {
                    siFactor = unitSystem.parse("Length/Time").getSIScaling();
                }

                for (size_t i = 0; i < water_vel_vals_.size(); ++i) {
                    water_vel_vals_[i] *= siFactor;
                }


                plyshlog_ref_conc_ = plyshlogTable.getRefPolymerConcentration();

                if (plyshlogTable.hasRefSalinity()) {
                    has_plyshlog_ref_salinity_ = true;
                    plyshlog_ref_salinity_ = plyshlogTable.getRefSalinity();
                } else {
                    has_plyshlog_ref_salinity_ = false;
                }

                if (plyshlogTable.hasRefTemperature()) {
                    has_plyshlog_ref_temp_ = true;
                    plyshlog_ref_temp_ = plyshlogTable.getRefTemperature();
                } else {
                    has_plyshlog_ref_temp_ = false;
                }
            }
        }

        double cMax() const;

        double mixParam() const;

        double rockDensity() const;

        double deadPoreVol() const;

        double resFactor() const;

        double cMaxAds() const;

        int adsIndex() const;

        /// indicate whehter PLYSHLOG is specified
        bool hasPlyshlog() const;

        /// the water velocity or water shear rate in PLYSHLOG table
        const std::vector<double>& shearWaterVelocity() const;

        /// the viscosity reduction factor PLYSHLOG table
        const std::vector<double>& shearViscosityReductionFactor() const;

        /// the reference polymer concentration in PLYSHLOG
        double plyshlogRefConc() const;

        /// indicate wheter reference salinity is specified in PLYSHLOG
        bool hasPlyshlogRefSalinity() const;

        /// indicate whether reference temperature is specified in PLYSHLOG
        bool hasPlyshlogRefTemp() const;

        /// the reference salinity in PLYSHLOG
        double plyshlogRefSalinity() const;

        /// the reference temperature in PLYSHLOG
        double plyshlogRefTemp() const;

        /// indicate whether SHRATE keyword is specified
        bool hasShrate() const;

        /// the value of SHRATE
        double shrate() const;

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

        void effectiveVisc(const double c, const double mu_w,
                                              double& mu_w_eff) const;

        void effectiveViscWithDer(const double c, const double visc
                                                     , double& mu_w_eff
                                                     , double dmu_w_eff_dc) const;

        void effectiveInvVisc(const double c, const double mu_w,
                                                 double& inv_mu_w_eff) const;

        void effectiveInvViscWithDer(const double c,
                                     const double mu_w,
                                     double& inv_mu_w_eff,
                                     double& dinv_mu_w_eff_dc) const;

        void effectiveInvPolyVisc(const double c,
                                  const double mu_w,
                                  double& inv_mu_p_eff) const;

        void effectiveInvPolyViscWithDer(const double c,
                                         const double mu_w,
                                         double& inv_mu_p_eff,
                                         double& d_inv_mu_p_eff_dc) const;

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

        /// Computing the shear multiplier based on the water velocity/shear rate with PLYSHLOG keyword
        bool computeShearMultLog(std::vector<double>& water_vel, std::vector<double>& visc_mult, std::vector<double>& shear_mult) const;

    private:
        double c_max_;
        double mix_param_;
        double rock_density_;
        double dead_pore_vol_;
        double res_factor_;
        double c_max_ads_;

        bool   has_plyshlog_;
        bool   has_shrate_;
        // Assuming NTPVT == 1 always due to the limitation of the parser
        // only one SHRATE value
        // TODO: to be extended later when parser is improved.
        double shrate_;
        AdsorptionBehaviour ads_index_;
        std::vector<double> c_vals_visc_;
        std::vector<double> visc_mult_vals_;
        std::vector<double> c_vals_ads_;
        std::vector<double> ads_vals_;
        std::vector<double> water_vel_vals_;
        std::vector<double> shear_vrf_vals_;

        double plyshlog_ref_conc_;
        double plyshlog_ref_salinity_;
        double plyshlog_ref_temp_;
        bool has_plyshlog_ref_salinity_;
        bool has_plyshlog_ref_temp_;


        void simpleAdsorptionBoth(double c, double& c_ads,
                                  double& dc_ads_dc, bool if_with_der) const;
        void adsorptionBoth(double c, double cmax,
                            double& c_ads, double& dc_ads_dc,
                            bool if_with_der) const;

        void effectiveInvViscBoth(const double c, const double mu_w,
                                  double& inv_mu_w_eff,
                                  double& dinv_mu_w_eff_dc, bool if_with_der) const;

        void effectiveInvPolyViscBoth(const double c,
                                      const double mu_w,
                                      double& inv_mu_p_eff,
                                      double& dinv_mu_p_eff_dc,
                                      const bool if_with_der) const;

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
