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
 * \copydoc Opm::SolventPvt
 */
#ifndef OPM_SOLVENT_PVT_HPP
#define OPM_SOLVENT_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvdsTable.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#endif

#include <vector>

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the "second"
 *        gas phase in the of ECL simulations with solvents.
 */
template <class Scalar>
class SolventPvt
{
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    explicit SolventPvt() = default;
    SolventPvt(const std::vector<Scalar>& solventReferenceDensity,
               const std::vector<TabulatedOneDFunction>& inverseSolventB,
               const std::vector<TabulatedOneDFunction>& solventMu,
               const std::vector<TabulatedOneDFunction>& inverseSolventBMu)
        : solventReferenceDensity_(solventReferenceDensity)
        , inverseSolventB_(inverseSolventB)
        , solventMu_(solventMu)
        , inverseSolventBMu_(inverseSolventBMu)
    {
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for "solvent gas" using an ECL deck.
     *
     * This method assumes that the deck features valid SDENSITY and PVDS keywords.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        const auto& pvdsTables = eclState.getTableManager().getPvdsTables();
        const auto& sdensityKeyword = deck.getKeyword("SDENSITY");

        assert(pvdsTables.size() == sdensityKeyword.size());

        size_t numRegions = pvdsTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefS = sdensityKeyword.getRecord(regionIdx).getItem("SOLVENT_DENSITY").getSIDouble(0);

            setReferenceDensity(regionIdx, rhoRefS);

            const auto& pvdsTable = pvdsTables.getTable<PvdsTable>(regionIdx);

            // say 99.97% of all time: "premature optimization is the root of all
            // evil". Eclipse does this "optimization" for apparently no good reason!
            std::vector<Scalar> invB(pvdsTable.numRows());
            const auto& Bg = pvdsTable.getFormationFactorColumn();
            for (unsigned i = 0; i < Bg.size(); ++ i) {
                invB[i] = 1.0/Bg[i];
            }

            size_t numSamples = invB.size();
            inverseSolventB_[regionIdx].setXYArrays(numSamples, pvdsTable.getPressureColumn(), invB);
            solventMu_[regionIdx].setXYArrays(numSamples, pvdsTable.getPressureColumn(), pvdsTable.getViscosityColumn());
        }

        initEnd();
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        solventReferenceDensity_.resize(numRegions);
        inverseSolventB_.resize(numRegions);
        inverseSolventBMu_.resize(numRegions);
        solventMu_.resize(numRegions);
    }

    /*!
     * \brief Initialize the reference density of the solvent gas for a given PVT region
     */
    void setReferenceDensity(unsigned regionIdx, Scalar rhoRefSolvent)
    { solventReferenceDensity_[regionIdx] = rhoRefSolvent; }

    /*!
     * \brief Initialize the viscosity of the solvent gas phase.
     *
     * This is a function of \f$(p_g)\f$...
     */
    void setSolventViscosity(unsigned regionIdx, const TabulatedOneDFunction& mug)
    { solventMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the function for the formation volume factor of solvent gas
     *
     * \param samplePoints A container of \f$(p_g, B_s)\f$ values
     */
    void setSolventFormationVolumeFactor(unsigned regionIdx, const SamplingPoints& samplePoints)
    {
        SamplingPoints tmp(samplePoints);
        auto it = tmp.begin();
        const auto& endIt = tmp.end();
        for (; it != endIt; ++ it)
            std::get<1>(*it) = 1.0/std::get<1>(*it);

        inverseSolventB_[regionIdx].setContainerOfTuples(tmp);
        assert(inverseSolventB_[regionIdx].monotonic());
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = solventMu_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the solvent
            // formation volume factor and its viscosity
            const auto& solventMu = solventMu_[regionIdx];
            const auto& invSolventB = inverseSolventB_[regionIdx];
            assert(solventMu.numSamples() == invSolventB.numSamples());

            std::vector<Scalar> pressureValues(solventMu.numSamples());
            std::vector<Scalar> invSolventBMuValues(solventMu.numSamples());
            for (unsigned pIdx = 0; pIdx < solventMu.numSamples(); ++pIdx) {
                pressureValues[pIdx] = invSolventB.xAt(pIdx);
                invSolventBMuValues[pIdx] = invSolventB.valueAt(pIdx) * (1.0/solventMu.valueAt(pIdx));
            }

            inverseSolventBMu_[regionIdx].setXYContainers(pressureValues, invSolventBMuValues);
        }
    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return solventReferenceDensity_.size(); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                                  const Evaluation& temperature OPM_UNUSED,
                                  const Evaluation& pressure) const
    {
        const Evaluation& invBg = inverseSolventB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseSolventBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Return the reference density the solvent phase for a given PVT region
     */
    Scalar referenceDensity(unsigned regionIdx) const
    { return solventReferenceDensity_[regionIdx]; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature OPM_UNUSED,
                                            const Evaluation& pressure) const
    { return inverseSolventB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    const std::vector<Scalar>& solventReferenceDensity() const
    { return solventReferenceDensity_; }

    const std::vector<TabulatedOneDFunction>& inverseSolventB() const
    { return inverseSolventB_; }

    const std::vector<TabulatedOneDFunction>& solventMu() const
    { return solventMu_; }

    const std::vector<TabulatedOneDFunction>& inverseSolventBMu() const
    { return inverseSolventBMu_; }

    bool operator==(const SolventPvt<Scalar>& data) const
    {
        return solventReferenceDensity_ == data.solventReferenceDensity_ &&
               inverseSolventB_ == data.inverseSolventB_ &&
               solventMu_ == data.solventMu_ &&
               inverseSolventBMu_ == data.inverseSolventBMu_;
    }

private:
    std::vector<Scalar> solventReferenceDensity_;
    std::vector<TabulatedOneDFunction> inverseSolventB_;
    std::vector<TabulatedOneDFunction> solventMu_;
    std::vector<TabulatedOneDFunction> inverseSolventBMu_;
};

} // namespace Opm

#endif
