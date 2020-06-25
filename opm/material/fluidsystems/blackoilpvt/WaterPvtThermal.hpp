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
 * \copydoc Opm::WaterPvtThermal
 */
#ifndef OPM_WATER_PVT_THERMAL_HPP
#define OPM_WATER_PVT_THERMAL_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {
template <class Scalar, bool enableThermal, bool enableBrine>
class WaterPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of water
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class WaterPvtThermal
{
public:
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef WaterPvtMultiplexer<Scalar, /*enableThermal=*/false, false> IsothermalPvt;

    WaterPvtThermal()
    {
        enableThermalDensity_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
        isothermalPvt_ = nullptr;
    }

    WaterPvtThermal(IsothermalPvt* isothermalPvt,
                    const std::vector<Scalar>& viscrefPress,
                    const std::vector<Scalar>& watdentRefTemp,
                    const std::vector<Scalar>& watdentCT1,
                    const std::vector<Scalar>& watdentCT2,
                    const std::vector<Scalar>& pvtwRefPress,
                    const std::vector<Scalar>& pvtwRefB,
                    const std::vector<Scalar>& pvtwCompressibility,
                    const std::vector<Scalar>& pvtwViscosity,
                    const std::vector<Scalar>& pvtwViscosibility,
                    const std::vector<TabulatedOneDFunction>& watvisctCurves,
                    const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                    bool enableThermalDensity,
                    bool enableThermalViscosity,
                    bool enableInternalEnergy)
        : isothermalPvt_(isothermalPvt)
        , viscrefPress_(viscrefPress)
        , watdentRefTemp_(watdentRefTemp)
        , watdentCT1_(watdentCT1)
        , watdentCT2_(watdentCT2)
        , pvtwRefPress_(pvtwRefPress)
        , pvtwRefB_(pvtwRefB)
        , pvtwCompressibility_(pvtwCompressibility)
        , pvtwViscosity_(pvtwViscosity)
        , pvtwViscosibility_(pvtwViscosibility)
        , watvisctCurves_(watvisctCurves)
        , internalEnergyCurves_(internalEnergyCurves)
        , enableThermalDensity_(enableThermalDensity)
        , enableThermalViscosity_(enableThermalViscosity)
        , enableInternalEnergy_(enableInternalEnergy)
    { }

    WaterPvtThermal(const WaterPvtThermal& data)
    { *this = data; }

    ~WaterPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the water PVT properties.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        //////
        // initialize the isothermal part
        //////
        isothermalPvt_ = new IsothermalPvt;
        isothermalPvt_->initFromState(eclState, schedule);

        //////
        // initialize the thermal part
        //////
        const auto& tables = eclState.getTableManager();

        enableThermalDensity_ = tables.WatDenT().size() > 0;
        enableThermalViscosity_ = tables.hasTables("WATVISCT");
        enableInternalEnergy_ = tables.hasTables("SPECHEAT");

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        if (enableThermalDensity_) {
            const auto& watDenT = tables.WatDenT();

            assert(watDenT.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& record = watDenT[regionIdx];

                watdentRefTemp_[regionIdx] = record.T0;
                watdentCT1_[regionIdx] = record.C1;
                watdentCT2_[regionIdx] = record.C2;
            }
        }

        if (enableThermalViscosity_) {
            if (tables.getViscrefTable().empty())
                throw std::runtime_error("VISCREF is required when WATVISCT is present");

            const auto& watvisctTables = tables.getWatvisctTables();
            const auto& viscrefTables = tables.getViscrefTable();

            const auto& pvtwTables = tables.getPvtwTable();

            assert(pvtwTables.size() == numRegions);
            assert(watvisctTables.size() == numRegions);
            assert(viscrefTables.size() == numRegions);

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                const auto& T = watvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& mu = watvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
                watvisctCurves_[regionIdx].setXYContainers(T, mu);

                viscrefPress_[regionIdx] = viscrefTables[regionIdx].reference_pressure;
            }

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                pvtwViscosity_[regionIdx] = pvtwTables[regionIdx].viscosity;
                pvtwViscosibility_[regionIdx] = pvtwTables[regionIdx].viscosibility;
            }
        }

        if (enableInternalEnergy_) {
            // the specific internal energy of liquid water. be aware that ecl only specifies the heat capacity
            // (via the SPECHEAT keyword) and we need to integrate it ourselfs to get the
            // internal energy
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& specHeatTable = tables.getSpecheatTables()[regionIdx];
                const auto& temperatureColumn = specHeatTable.getColumn("TEMPERATURE");
                const auto& cvWaterColumn = specHeatTable.getColumn("CV_WATER");

                std::vector<double> uSamples(temperatureColumn.size());

                Scalar u = temperatureColumn[0]*cvWaterColumn[0];
                for (size_t i = 0;; ++i) {
                    uSamples[i] = u;

                    if (i >= temperatureColumn.size() - 1)
                        break;

                    // integrate to the heat capacity from the current sampling point to the next
                    // one. this leads to a quadratic polynomial.
                    Scalar c_v0 = cvWaterColumn[i];
                    Scalar c_v1 = cvWaterColumn[i + 1];
                    Scalar T0 = temperatureColumn[i];
                    Scalar T1 = temperatureColumn[i + 1];
                    u += 0.5*(c_v0 + c_v1)*(T1 - T0);
                }

                internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
            }
        }
    }
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    {
        pvtwRefPress_.resize(numRegions);
        pvtwRefB_.resize(numRegions);
        pvtwCompressibility_.resize(numRegions);
        pvtwViscosity_.resize(numRegions);
        pvtwViscosibility_.resize(numRegions);
        viscrefPress_.resize(numRegions);
        watvisctCurves_.resize(numRegions);
        watdentRefTemp_.resize(numRegions);
        watdentCT1_.resize(numRegions);
        watdentCT2_.resize(numRegions);
        internalEnergyCurves_.resize(numRegions);
    }

    /*!
     * \brief Finish initializing the thermal part of the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the water phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff the viscosity of the water phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return pvtwRefPress_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of water given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure OPM_UNUSED) const
    {
        if (!enableInternalEnergy_)
            throw std::runtime_error("Requested the internal energy of oil but it is disabled");

        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& saltconcentration) const
    {
        const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, saltconcentration);
        if (!enableThermalViscosity())
            return isothermalMu;

        Scalar x = -pvtwViscosibility_[regionIdx]*(viscrefPress_[regionIdx] - pvtwRefPress_[regionIdx]);
        Scalar muRef = pvtwViscosity_[regionIdx]/(1.0 + x + 0.5*x*x);

        // compute the viscosity deviation due to temperature
        const auto& muWatvisct = watvisctCurves_[regionIdx].eval(temperature, true);
        return isothermalMu * muWatvisct/muRef;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& saltconcentration) const
    {
        if (!enableThermalDensity())
            return isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, saltconcentration);

        Scalar BwRef = pvtwRefB_[regionIdx];
        Scalar TRef = watdentRefTemp_[regionIdx];
        const Evaluation& X = pvtwCompressibility_[regionIdx]*(pressure - pvtwRefPress_[regionIdx]);
        Scalar cT1 = watdentCT1_[regionIdx];
        Scalar cT2 = watdentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        // this is inconsistent with the density calculation of water in the isothermal
        // case (it misses the quadratic pressure term), but it is the equation given in
        // the documentation.
        return 1.0/(((1 - X)*(1 + cT1*Y + cT2*Y*Y))*BwRef);
    }

    const IsothermalPvt* isoThermalPvt() const
    { return isothermalPvt_; }

    const std::vector<Scalar>& viscrefPress() const
    { return viscrefPress_; }

    const std::vector<Scalar>& watdentRefTemp() const
    { return watdentRefTemp_; }

    const std::vector<Scalar>& watdentCT1() const
    { return watdentCT1_; }

    const std::vector<Scalar>& watdentCT2() const
    { return watdentCT2_; }

    const std::vector<Scalar>& pvtwRefPress() const
    { return pvtwRefPress_; }

    const std::vector<Scalar>& pvtwRefB() const
    { return pvtwRefB_; }

    const std::vector<Scalar>& pvtwCompressibility() const
    { return pvtwCompressibility_; }

    const std::vector<Scalar>& pvtwViscosity() const
    { return pvtwViscosity_; }

    const std::vector<Scalar>& pvtwViscosibility() const
    { return pvtwViscosibility_; }

    const std::vector<TabulatedOneDFunction>& watvisctCurves() const
    { return watvisctCurves_; }

    const std::vector<TabulatedOneDFunction> internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    bool operator==(const WaterPvtThermal<Scalar>& data) const
    {
        if (isothermalPvt_ && !data.isothermalPvt_)
            return false;
        if (!isothermalPvt_ && data.isothermalPvt_)
            return false;

        return (!this->isoThermalPvt() ||
               (*this->isoThermalPvt() == *data.isoThermalPvt())) &&
               this->viscrefPress() == data.viscrefPress() &&
               this->watdentRefTemp() == data.watdentRefTemp() &&
               this->watdentCT1() == data.watdentCT1() &&
               this->watdentCT2() == data.watdentCT2() &&
               this->pvtwRefPress() == data.pvtwRefPress() &&
               this->pvtwRefB() == data.pvtwRefB() &&
               this->pvtwCompressibility() == data.pvtwCompressibility() &&
               this->pvtwViscosity() == data.pvtwViscosity() &&
               this->pvtwViscosibility() == data.pvtwViscosibility() &&
               this->watvisctCurves() == data.watvisctCurves() &&
               this->internalEnergyCurves() == data.internalEnergyCurves() &&
               this->enableThermalDensity() == data.enableThermalDensity() &&
               this->enableThermalViscosity() == data.enableThermalViscosity() &&
               this->enableInternalEnergy() == data.enableInternalEnergy();
    }

    WaterPvtThermal<Scalar>& operator=(const WaterPvtThermal<Scalar>& data)
    {
        if (data.isothermalPvt_)
            isothermalPvt_ = new IsothermalPvt(*data.isothermalPvt_);
        else
            isothermalPvt_ = nullptr;
        viscrefPress_ = data.viscrefPress_;
        watdentRefTemp_ = data.watdentRefTemp_;
        watdentCT1_ = data.watdentCT1_;
        watdentCT2_ = data.watdentCT2_;
        pvtwRefPress_ = data.pvtwRefPress_;
        pvtwRefB_ = data.pvtwRefB_;
        pvtwCompressibility_ = data.pvtwCompressibility_;
        pvtwViscosity_ = data.pvtwViscosity_;
        pvtwViscosibility_ = data.pvtwViscosibility_;
        watvisctCurves_ = data.watvisctCurves_;
        internalEnergyCurves_ = data.internalEnergyCurves_;
        enableThermalDensity_ = data.enableThermalDensity_;
        enableThermalViscosity_ = data.enableThermalViscosity_;
        enableInternalEnergy_ = data.enableInternalEnergy_;

        return *this;
    }

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence. We need to store one
    // value per PVT region.
    std::vector<Scalar> viscrefPress_;

    std::vector<Scalar> watdentRefTemp_;
    std::vector<Scalar> watdentCT1_;
    std::vector<Scalar> watdentCT2_;

    std::vector<Scalar> pvtwRefPress_;
    std::vector<Scalar> pvtwRefB_;
    std::vector<Scalar> pvtwCompressibility_;
    std::vector<Scalar> pvtwViscosity_;
    std::vector<Scalar> pvtwViscosibility_;

    std::vector<TabulatedOneDFunction> watvisctCurves_;

    // piecewise linear curve representing the internal energy of water
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Opm

#endif
