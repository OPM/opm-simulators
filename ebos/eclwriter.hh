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
 * \copydoc Ewoms::EclWriter
 */
#ifndef EWOMS_ECL_WRITER_HH
#define EWOMS_ECL_WRITER_HH

#include <opm/material/densead/Evaluation.hpp>

#include "collecttoiorank.hh"
#include "ecloutputblackoilmodule.hh"

#include <ewoms/disc/ecfv/ecfvdiscretization.hh>
#include <ewoms/io/baseoutputwriter.hh>
#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/common/Valgrind.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <boost/algorithm/string.hpp>

#include <list>
#include <utility>
#include <string>
#include <limits>
#include <sstream>
#include <fstream>
#include <type_traits>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(EnableEclOutput);
}

template <class TypeTag>
class EclWriter;


/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to opm-output.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you will have to
 *   have the OPM module opm-output.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "opm-grid". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an ECL output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
template <class TypeTag>
class EclWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) GridManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;


    typedef CollectDataToIORank< GridManager > CollectDataToIORankType;

    typedef std::vector<Scalar> ScalarBuffer;


public:
    EclWriter(const Simulator& simulator)
        : simulator_(simulator)
        , eclOutputModule_(simulator)
        , collectToIORank_( simulator_.gridManager() )
    {
        Grid globalGrid = simulator_.gridManager().grid();
        globalGrid.switchToGlobalView();
        eclIO_.reset(new Opm::EclipseIO(simulator_.gridManager().eclState(),
                                        Opm::UgGridHelpers::createEclipseGrid( globalGrid , simulator_.gridManager().eclState().getInputGrid() ),
                                        simulator_.gridManager().schedule(),
                                        simulator_.gridManager().summaryConfig()));
    }

    ~EclWriter()
    { }

    void setEclIO(std::unique_ptr<Opm::EclipseIO>&& eclIO) {
        eclIO_ = std::move(eclIO);
    }

    const Opm::EclipseIO& eclIO() const
    {return *eclIO_;}

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */
    void writeOutput(const Opm::data::Wells& dw, Scalar t, bool substep, Scalar totalSolverTime, Scalar nextstep, const Opm::data::Solution& fip)
    {

        #if !HAVE_OPM_OUTPUT
                OPM_THROW(std::runtime_error,
                          "Opm-output must be available to write ECL output!");
        #else

        int episodeIdx = simulator_.episodeIndex() + 1;
        const auto& gridView = simulator_.gridManager().gridView();
        int numElements = gridView.size(/*codim=*/0);
        bool log = collectToIORank_.isIORank();
        eclOutputModule_.allocBuffers(numElements, episodeIdx, simulator_.gridManager().eclState().getRestartConfig(), log);

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView.template begin</*codim=*/0>();
        const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            eclOutputModule_.processElement(elemCtx);
        }
        eclOutputModule_.outputErrorLog();

        // collect all data to I/O rank and assign to sol
        Opm::data::Solution localCellData = fip;
        eclOutputModule_.assignToSolution(localCellData);
        if (collectToIORank_.isParallel())
            collectToIORank_.collect(localCellData);

        // write output on I/O rank
        if (collectToIORank_.isIORank()) {

            std::map<std::string, std::vector<double>> extraRestartData;
            std::map<std::string, double> miscSummaryData;

            // Add suggested next timestep to extra data.
            extraRestartData["OPMEXTRA"] = std::vector<double>(1, nextstep);

            // Add TCPU if simulatorReport is not defaulted.
            if (totalSolverTime != 0.0) {
                miscSummaryData["TCPU"] = totalSolverTime;
            }

            const Opm::data::Solution& cellData = collectToIORank_.isParallel() ? collectToIORank_.globalCellData() : localCellData;
            eclIO_->writeTimeStep(episodeIdx,
                                  substep,
                                  t,
                                  cellData,
                                  dw,
                                  miscSummaryData,
                                  extraRestartData,
                                  false);
        }

#endif
    }

    void restartBegin() {
        std::map<std::string, Opm::RestartKey> solution_keys {{"PRESSURE" , Opm::RestartKey(Opm::UnitSystem::measure::pressure)},
                                                         {"SWAT" , Opm::RestartKey(Opm::UnitSystem::measure::identity)},
                                                         {"SGAS" , Opm::RestartKey(Opm::UnitSystem::measure::identity)},
                                                         {"TEMP" , Opm::RestartKey(Opm::UnitSystem::measure::temperature)},
                                                         {"RS" , Opm::RestartKey(Opm::UnitSystem::measure::gas_oil_ratio)},
                                                         {"RV" , Opm::RestartKey(Opm::UnitSystem::measure::oil_gas_ratio)},
                                                         {"SOMAX", {Opm::UnitSystem::measure::identity, false}},
                                                         {"PCSWM_OW", {Opm::UnitSystem::measure::identity, false}},
                                                         {"KRNSW_OW", {Opm::UnitSystem::measure::identity, false}},
                                                         {"PCSWM_GO", {Opm::UnitSystem::measure::identity, false}},
                                                         {"KRNSW_GO", {Opm::UnitSystem::measure::identity, false}}};

        std::map<std::string, bool> extra_keys {
            {"OPMEXTRA" , false}
        };

        unsigned episodeIdx = simulator_.episodeIndex();
        const auto& gridView = simulator_.gridManager().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        eclOutputModule_.allocBuffers(numElements, episodeIdx, simulator_.gridManager().eclState().getRestartConfig(), false);

        auto restart_values = eclIO_->loadRestart(solution_keys, extra_keys);
        for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
            unsigned globalIdx = collectToIORank_.localIdxToGlobalIdx(elemIdx);
            eclOutputModule_.setRestart(restart_values.solution, elemIdx, globalIdx);
        }
    }


    const EclOutputBlackOilModule<TypeTag>& eclOutputModule() const {
        return eclOutputModule_;
    }


private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    const Simulator& simulator_;
    EclOutputBlackOilModule<TypeTag> eclOutputModule_;
    CollectDataToIORankType collectToIORank_;
    std::unique_ptr<Opm::EclipseIO> eclIO_;

};
} // namespace Ewoms

#endif
