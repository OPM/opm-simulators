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
 * \copydoc Opm::EclProblem
 */
#ifndef EWOMS_GENERIC_ECL_PROBLEM_HH
#define EWOMS_GENERIC_ECL_PROBLEM_HH

#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#include <array>
#include <string>
#include <vector>

namespace Opm {

class Deck;
class EclipseState;
class Schedule;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template<class GridView, class FluidSystem, class Scalar>
class EclGenericProblem
{
public:
    using TabulatedTwoDFunction = UniformXTabulated2DFunction<Scalar>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    struct RockParams {
        Scalar referencePressure;
        Scalar compressibility;
    };

    EclGenericProblem(const EclipseState& eclState,
                      const Schedule& schedule,
                      const GridView& gridView);

    /*!
     * \copydoc FvBaseProblem::helpPreamble
     */
    static std::string helpPreamble(int,
                                    const char **argv);

    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription();

    /*!
     * \brief Specifies the string returned by briefDescription()
     *
     * This string appears in the usage message.
     */
    static void setBriefDescription(const std::string& msg)
    { briefDescription_ = msg; }

    /*!
     * \brief Returns an element's historic maximum water phase saturation that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is used for output of the maximum water saturation used as input
     * for water induced rock compation ROCK2D/ROCK2DTR.
     */
    Scalar maxWaterSaturation(unsigned globalDofIdx) const;

    /*!
     * \brief Returns an element's historic minimum pressure of the oil phase that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is used for output of the minimum pressure used as input
     * for the irreversible rock compation option.
     */
    Scalar minOilPressure(unsigned globalDofIdx) const;

    /*!
     * \brief Get the pressure of the overburden.
     *
     * This method is mainly for output.
     */
    Scalar overburdenPressure(unsigned elementIdx) const;

    /*!
     * \brief Returns the porosity of an element
     *
     * The reference porosity of an element is the porosity of the medium before modified
     * by the current solution. Note that this method is *not* part of the generic eWoms
     * problem API because it would bake the assumption that only the elements are the
     * degrees of freedom into the interface.
     */
    Scalar referencePorosity(unsigned elementIdx, unsigned timeIdx) const
    { return referencePorosity_[timeIdx][elementIdx]; }

    /*!
     * \brief Sets the porosity of an element
     *
     */
    void setPorosity(Scalar poro, unsigned elementIdx, unsigned timeIdx = 0)
    { referencePorosity_[timeIdx][elementIdx] = poro; }

    /*!
     * \brief Returns the initial solvent saturation for a given a cell index
     */
    Scalar solventSaturation(unsigned elemIdx) const;

    /*!
     * \brief Returns the dynamic drsdt convective mixing value
     */
    Scalar drsdtcon(unsigned elemIdx, int episodeIdx) const;

    /*!
     * \brief Returns the initial polymer concentration for a given a cell index
     */
    Scalar  polymerConcentration(unsigned elemIdx) const;

    /*!
    * \brief Returns the polymer molecule weight for a given cell index
    */
    // TODO: remove this function if not called
    Scalar polymerMolecularWeight(const unsigned elemIdx) const;

    /*!
     * \brief Returns the initial microbial concentration for a given a cell index
     */
    Scalar  microbialConcentration(unsigned elemIdx) const;

    /*!
     * \brief Returns the initial oxygen concentration for a given a cell index
     */
    Scalar  oxygenConcentration(unsigned elemIdx) const;

    /*!
     * \brief Returns the initial urea concentration for a given a cell index
     */
    Scalar  ureaConcentration(unsigned elemIdx) const;

    /*!
     * \brief Returns the initial biofilm concentration for a given a cell index
     */
    Scalar  biofilmConcentration(unsigned elemIdx) const;

    /*!
     * \brief Returns the initial calcite concentration for a given a cell index
     */
    Scalar  calciteConcentration(unsigned elemIdx) const;

    /*!
     * \brief Returns the index the relevant PVT region given a cell index
     */
    unsigned pvtRegionIndex(unsigned elemIdx) const;

//    const std::vector<int>& pvtRegionArray() const
//    { return pvtnum_; }

    /*!
     * \brief Returns the index the relevant saturation function region given a cell index
     */
    unsigned satnumRegionIndex(unsigned elemIdx) const;

    /*!
     * \brief Returns the index the relevant MISC region given a cell index
     */
    unsigned miscnumRegionIndex(unsigned elemIdx) const;

    /*!
     * \brief Returns the index the relevant PLMIXNUM (for polymer module) region given a cell index
     */
    unsigned plmixnumRegionIndex(unsigned elemIdx) const;

    /*!
     * \brief Returns the max polymer adsorption value
     */
    Scalar maxPolymerAdsorption(unsigned elemIdx) const;

    /*!
     * Direct access to rock compressibility.
     *
     * While the above overload could be implemented in terms of this method,
     * that would require always looking up the global space index, which
     * is not always needed.
     */
    Scalar rockCompressibility(unsigned globalSpaceIdx) const;

    /*!
     * Direct access to rock reference pressure.
     *
     * While the above overload could be implemented in terms of this method,
     * that would require always looking up the global space index, which
     * is not always needed.
     */
    Scalar rockReferencePressure(unsigned globalSpaceIdx) const;

    /*!
     * \brief Direct indexed access to the porosity.
     *
     * For the EclProblem, this method is identical to referencePorosity(). The intensive
     * quantities object may apply various multipliers (e.g. ones which model rock
     * compressibility and water induced rock compaction) to it which depend on the
     * current physical conditions.
     */
    Scalar porosity(unsigned globalSpaceIdx, unsigned timeIdx) const;

    /*!
     * \brief Returns the minimum allowable size of a time step.
     */
    Scalar minTimeStepSize() const
    { return minTimeStepSize_; }

    /*!
     * \brief Returns the maximum number of subsequent failures for the time integration
     *        before giving up.
     */
    unsigned maxTimeIntegrationFailures() const
    { return maxFails_; }

    bool vapparsActive(int episodeIdx) const;

protected:
    bool drsdtActive_(int episodeIdx) const;
    bool drvdtActive_(int episodeIdx) const;
    bool drsdtConvective_(int episodeIdx) const;

    void initFluidSystem_();
    void initDRSDT_(size_t numDof,
                    int episodeIdx);

    /*!
     * \brief Always returns true. The ecl output writer takes care of the rest
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Returns true if an eWoms restart file should be written to disk.
     *
     * The EclProblem does not write any restart files using the ad-hoc format, only ones
     * using the ECL format.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    bool beginEpisode_(bool enableExperiments,
                       int episodeIdx);
    void beginTimeStep_(bool enableExperiments,
                        int episodeIdx,
                        int timeStepIndex,
                        Scalar startTime,
                        Scalar time,
                        Scalar timeStepSize,
                        Scalar endTime);

    void readRockParameters_(const std::vector<Scalar>& cellCenterDepths);
    void readRockCompactionParameters_();

    void readBlackoilExtentionsInitialConditions_(size_t numDof,
                                                  bool enableSolvent,
                                                  bool enablePolymer,
                                                  bool enablePolymerMolarWeight,
                                                  bool enableMICP);

    void updatePvtnum_();
    void updateSatnum_();
    void updateMiscnum_();
    void updatePlmixnum_();
    void updateKrnum_();

    const EclipseState& eclState_;
    const Schedule& schedule_;
    const GridView& gridView_;

    static inline std::string briefDescription_;
    std::array<std::vector<Scalar>, 2> referencePorosity_;

    std::vector<int> pvtnum_;
    std::vector<unsigned short> satnum_;
    std::vector<unsigned short> miscnum_;
    std::vector<unsigned short> plmixnum_;
    std::vector<unsigned short> krnumx_;
    std::vector<unsigned short> krnumy_;
    std::vector<unsigned short> krnumz_;

    std::vector<RockParams> rockParams_;
    std::vector<unsigned short> rockTableIdx_;
    std::vector<TabulatedTwoDFunction> rockCompPoroMultWc_;
    std::vector<TabulatedTwoDFunction> rockCompTransMultWc_;
    std::vector<TabulatedFunction> rockCompPoroMult_;
    std::vector<TabulatedFunction> rockCompTransMult_;

    std::vector<Scalar> maxOilSaturation_;
    std::vector<Scalar> maxPolymerAdsorption_;
    std::vector<Scalar> maxWaterSaturation_;
    std::vector<Scalar> minOilPressure_;
    std::vector<Scalar> overburdenPressure_;
    std::vector<Scalar> polymerConcentration_;
    std::vector<Scalar> polymerMoleWeight_; // polymer molecular weight
    std::vector<Scalar> solventSaturation_;
    std::vector<Scalar> microbialConcentration_;
    std::vector<Scalar> oxygenConcentration_;
    std::vector<Scalar> ureaConcentration_;
    std::vector<Scalar> biofilmConcentration_;
    std::vector<Scalar> calciteConcentration_;

    std::vector<Scalar> lastRv_;
    std::vector<Scalar> maxDRv_;

    std::vector<Scalar> convectiveDrs_;
    std::vector<Scalar> lastRs_;
    std::vector<Scalar> maxDRs_;
    std::vector<bool> dRsDtOnlyFreeGas_; // apply the DRSDT rate limit only to cells that exhibit free gas

    // time stepping parameters
    bool enableTuning_;
    Scalar initialTimeStepSize_;
    Scalar maxTimeStepAfterWellEvent_;
    Scalar maxTimeStepSize_;
    Scalar restartShrinkFactor_;
    unsigned maxFails_;
    Scalar minTimeStepSize_;

private:
    template<class T>
    void updateNum(const std::string& name, std::vector<T>& numbers);
};

} // namespace Opm

#endif
