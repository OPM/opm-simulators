/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#ifndef OPM_INITSTATEEQUIL_HEADER_INCLUDED
#define OPM_INITSTATEEQUIL_HEADER_INCLUDED

#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/simulator/EquilibrationHelpers.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/RegionMapping.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RsvdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RvvdTable.hpp>

#include <array>
#include <cassert>
#include <utility>
#include <vector>

/**
 * \file
 * Facilities for an ECLIPSE-style equilibration-based
 * initialisation scheme (keyword 'EQUIL').
 */
struct UnstructuredGrid;

namespace Opm
{

    /**
     * Compute initial state by an equilibration procedure.
     *
     * The following state fields are modified:
     *   pressure(),
     *   saturation(),
     *   surfacevol(),
     *   gasoilratio(),
     *   rv().
     *
     * \param[in] grid     Grid.
     * \param[in] props    Property object, pvt and capillary properties are used.
     * \param[in] deck     Simulation deck, used to obtain EQUIL and related data.
     * \param[in] gravity  Acceleration of gravity, assumed to be in Z direction.
     */
    template<class Grid>
    void initStateEquil(const Grid& grid,
                        const BlackoilPropertiesInterface& props,
                        const Opm::DeckConstPtr deck,
                        const Opm::EclipseStateConstPtr eclipseState,
                        const double gravity,
                        BlackoilState& state);


    /**
     * Types and routines that collectively implement a basic
     * ECLIPSE-style equilibration-based initialisation scheme.
     *
     * This namespace is intentionally nested to avoid name clashes
     * with other parts of OPM.
     */
    namespace EQUIL {

        /**
         * Compute initial phase pressures by means of equilibration.
         *
         * This function uses the information contained in an
         * equilibration record (i.e., depths and pressurs) as well as
         * a density calculator and related data to vertically
         * integrate the phase pressure ODE
         * \f[
         * \frac{\mathrm{d}p_{\alpha}}{\mathrm{d}z} =
         * \rho_{\alpha}(z,p_{\alpha})\cdot g
         * \f]
         * in which \f$\rho_{\alpha}$ denotes the fluid density of
         * fluid phase \f$\alpha\f$, \f$p_{\alpha}\f$ is the
         * corresponding phase pressure, \f$z\f$ is the depth and
         * \f$g\f$ is the acceleration due to gravity (assumed
         * directed downwords, in the positive \f$z\f$ direction).
         *
         * \tparam Region Type of an equilibration region information
         *                base.  Typically an instance of the EquilReg
         *                class template.
         *
         * \tparam CellRange Type of cell range that demarcates the
         *                cells pertaining to the current
         *                equilibration region.  Must implement
         *                methods begin() and end() to bound the range
         *                as well as provide an inner type,
         *                const_iterator, to traverse the range.
         *
         * \param[in] G     Grid.
         * \param[in] reg   Current equilibration region.
         * \param[in] cells Range that spans the cells of the current
         *                  equilibration region.
         * \param[in] grav  Acceleration of gravity.
         *
         * \return Phase pressures, one vector for each active phase,
         * of pressure values in each cell in the current
         * equilibration region.
         */
        template <class Grid, class Region, class CellRange>
        std::vector< std::vector<double> >
        phasePressures(const Grid&             G,
                       const Region&           reg,
                       const CellRange&        cells,
                       const double            grav = unit::gravity);



        /**
         * Compute initial phase saturations by means of equilibration.
         *
         * \tparam Region Type of an equilibration region information
         *                base.  Typically an instance of the EquilReg
         *                class template.
         *
         * \tparam CellRange Type of cell range that demarcates the
         *                cells pertaining to the current
         *                equilibration region.  Must implement
         *                methods begin() and end() to bound the range
         *                as well as provide an inner type,
         *                const_iterator, to traverse the range.
         *
         * \param[in] reg             Current equilibration region.
         * \param[in] cells           Range that spans the cells of the current
         *                            equilibration region.
         * \param[in] props           Property object, needed for capillary functions.
         * \param[in] phase_pressures Phase pressures, one vector for each active phase,
         *                            of pressure values in each cell in the current
         *                            equilibration region.
         * \return                    Phase saturations, one vector for each phase, each containing
         *                            one saturation value per cell in the region.
         */
        template <class Grid, class Region, class CellRange>
        std::vector< std::vector<double> >
        phaseSaturations(const Grid&             grid,
                         const Region&           reg,
                         const CellRange&        cells,
                         BlackoilPropertiesInterface& props,
                         const std::vector<double> swat_init,
                         std::vector< std::vector<double> >& phase_pressures);



        /**
         * Compute initial Rs values.
         *
         * \tparam CellRangeType Type of cell range that demarcates the
         *                cells pertaining to the current
         *                equilibration region.  Must implement
         *                methods begin() and end() to bound the range
         *                as well as provide an inner type,
         *                const_iterator, to traverse the range.
         *
         * \param[in] grid            Grid.
         * \param[in] cells           Range that spans the cells of the current
         *                            equilibration region.
         * \param[in] oil_pressure    Oil pressure for each cell in range.
         * \param[in] temperature     Temperature for each cell in range.
         * \param[in] rs_func         Rs as function of pressure and depth.
         * \return                    Rs values, one for each cell in the 'cells' range.
         */
        template <class Grid, class CellRangeType>
        std::vector<double> computeRs(const Grid& grid,
                                      const CellRangeType& cells,
                                      const std::vector<double> oil_pressure,
                                      const std::vector<double>& temperature,
                                      const Miscibility::RsFunction& rs_func,
                                      const std::vector<double> gas_saturation);

        namespace DeckDependent {
            inline
            std::vector<EquilRecord>
            getEquil(const Opm::EclipseState& state)
            {
                const auto& init = *state.getInitConfig();

                if( !init.hasEquil() ) {
                    OPM_THROW(std::domain_error, "Deck does not provide equilibration data.");
                }

                const auto& equil = init.getEquil();
                return { equil.begin(), equil.end() };
            }

            template<class Grid>
            inline
            std::vector<int>
            equilnum(const Opm::DeckConstPtr deck,
                     const Opm::EclipseStateConstPtr eclipseState,
                     const Grid&  G   )
            {
                std::vector<int> eqlnum;
                if (deck->hasKeyword("EQLNUM")) {
                    const int nc = UgGridHelpers::numCells(G);
                    eqlnum.resize(nc);
                    const std::vector<int>& e = 
                        eclipseState->getIntGridProperty("EQLNUM")->getData();                    
                    const int* gc = UgGridHelpers::globalCell(G);
                    for (int cell = 0; cell < nc; ++cell) {
                        const int deck_pos = (gc == NULL) ? cell : gc[cell];
                        eqlnum[cell] = e[deck_pos] - 1;
                    }
                }
                else {
                    // No explicit equilibration region.
                    // All cells in region zero.
                    eqlnum.assign(UgGridHelpers::numCells(G), 0);
                }

                return eqlnum;
            }


            class InitialStateComputer {
            public:
                template<class Grid>
                InitialStateComputer(BlackoilPropertiesInterface& props,
                                     const Opm::DeckConstPtr            deck,
                                     const Opm::EclipseStateConstPtr eclipseState,
                                     const Grid&                        G    ,
                                     const double                       grav = unit::gravity)
                    : pp_(props.numPhases(),
                          std::vector<double>(UgGridHelpers::numCells(G))),
                      sat_(props.numPhases(),
                          std::vector<double>(UgGridHelpers::numCells(G))),
                      rs_(UgGridHelpers::numCells(G)),
                      rv_(UgGridHelpers::numCells(G))
                {
                    // Get the equilibration records.
                    const std::vector<EquilRecord> rec = getEquil(*eclipseState);
                    std::shared_ptr<const TableManager> tables = eclipseState->getTableManager();
                    // Create (inverse) region mapping.
                    const RegionMapping<> eqlmap(equilnum(deck, eclipseState, G)); 

                    // Create Rs functions.
                    rs_func_.reserve(rec.size());
                    if (deck->hasKeyword("DISGAS")) {                    
                        const TableContainer& rsvdTables = tables->getRsvdTables();
                        for (size_t i = 0; i < rec.size(); ++i) {
                            const int cell = *(eqlmap.cells(i).begin());                   
                            if (!rec[i].liveOilInitConstantRs()) {
                                if (rsvdTables.size() <= 0 ) {
                                    OPM_THROW(std::runtime_error, "Cannot initialise: RSVD table not available.");
                                }
                                const RsvdTable& rsvdTable = rsvdTables.getTable<RsvdTable>(i);
                                std::vector<double> depthColumn = rsvdTable.getColumn("DEPTH").vectorCopy();
                                std::vector<double> rsColumn = rsvdTable.getColumn("RS").vectorCopy();
                                rs_func_.push_back(std::make_shared<Miscibility::RsVD>(props,
                                                                                        cell,
                                                                                        depthColumn , rsColumn));
                            } else {
                                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                                    OPM_THROW(std::runtime_error,
                                              "Cannot initialise: when no explicit RSVD table is given, \n"
                                              "datum depth must be at the gas-oil-contact. "
                                              "In EQUIL region " << (i + 1) << "  (counting from 1), this does not hold.");
                                }
                                const double p_contact = rec[i].datumDepthPressure();
                                const double T_contact = 273.15 + 20; // standard temperature for now
                                rs_func_.push_back(std::make_shared<Miscibility::RsSatAtContact>(props, cell, p_contact, T_contact));
                            }
                        }
                    } else {
                        for (size_t i = 0; i < rec.size(); ++i) {
                            rs_func_.push_back(std::make_shared<Miscibility::NoMixing>());
                        }
                    }                    

                    rv_func_.reserve(rec.size());
                    if (deck->hasKeyword("VAPOIL")) {                    
                        const TableContainer& rvvdTables = tables->getRvvdTables();
                        for (size_t i = 0; i < rec.size(); ++i) {
                            const int cell = *(eqlmap.cells(i).begin());                   
                            if (!rec[i].wetGasInitConstantRv()) {
                                if (rvvdTables.size() <= 0) { 
                                    OPM_THROW(std::runtime_error, "Cannot initialise: RVVD table not available.");
                                }

                                const RvvdTable& rvvdTable = rvvdTables.getTable<RvvdTable>(i);
                                std::vector<double> depthColumn = rvvdTable.getColumn("DEPTH").vectorCopy();
                                std::vector<double> rvColumn = rvvdTable.getColumn("RV").vectorCopy();

                                rv_func_.push_back(std::make_shared<Miscibility::RvVD>(props,
                                                                                        cell,
                                                                                        depthColumn , rvColumn));

                            } else {
                                if (rec[i].gasOilContactDepth() != rec[i].datumDepth()) {
                                    OPM_THROW(std::runtime_error,
                                              "Cannot initialise: when no explicit RVVD table is given, \n"
                                              "datum depth must be at the gas-oil-contact. "
                                              "In EQUIL region " << (i + 1) << "  (counting from 1), this does not hold.");
                                }
                                const double p_contact = rec[i].datumDepthPressure() + rec[i].gasOilContactCapillaryPressure();
                                const double T_contact = 273.15 + 20; // standard temperature for now
                                rv_func_.push_back(std::make_shared<Miscibility::RvSatAtContact>(props, cell, p_contact, T_contact));
                            }
                        }
                    } else {
                        for (size_t i = 0; i < rec.size(); ++i) {
                            rv_func_.push_back(std::make_shared<Miscibility::NoMixing>());
                        }
                    }
                    

                    // Check for presence of kw SWATINIT
                    if (deck->hasKeyword("SWATINIT")) {
                        const std::vector<double>& swat_init = eclipseState->getDoubleGridProperty("SWATINIT")->getData();
                        const int nc = UgGridHelpers::numCells(G);
                        swat_init_.resize(nc);
                        const int* gc = UgGridHelpers::globalCell(G);
                        for (int c = 0; c < nc; ++c) {
                            const int deck_pos = (gc == NULL) ? c : gc[c];
                            swat_init_[c] = swat_init[deck_pos];
                        }
                    }

                    // Compute pressures, saturations, rs and rv factors.
                    calcPressSatRsRv(eqlmap, rec, props, G, grav);

                    // Modify oil pressure in no-oil regions so that the pressures of present phases can
                    // be recovered from the oil pressure and capillary relations.
                }

                typedef std::vector<double> Vec;
                typedef std::vector<Vec>    PVec; // One per phase.

                const PVec& press() const { return pp_; }
                const PVec& saturation() const { return sat_; }
                const Vec& rs() const { return rs_; }
                const Vec& rv() const { return rv_; }

            private:
                typedef DensityCalculator<BlackoilPropertiesInterface> RhoCalc;
                typedef EquilReg<RhoCalc> EqReg;

                std::vector< std::shared_ptr<Miscibility::RsFunction> > rs_func_;
                std::vector< std::shared_ptr<Miscibility::RsFunction> > rv_func_;

                PVec pp_;
                PVec sat_;
                Vec rs_;
                Vec rv_;
                Vec swat_init_;

                template <class RMap, class Grid>
                void
                calcPressSatRsRv(const RMap&                       reg  ,
                                 const std::vector< EquilRecord >& rec  ,
                                 Opm::BlackoilPropertiesInterface& props,
                                 const Grid&                       G    ,
                                 const double grav)
                {
                    for (const auto& r : reg.activeRegions()) {
                        const auto& cells = reg.cells(r);
                        const int repcell = *cells.begin();

                        const RhoCalc calc(props, repcell);
                        const EqReg eqreg(rec[r], calc,
                                          rs_func_[r], rv_func_[r],
                                          props.phaseUsage());
                   
                        PVec pressures = phasePressures(G, eqreg, cells, grav);
                        const std::vector<double>& temp = temperature(G, eqreg, cells);

                        const PVec sat = phaseSaturations(G, eqreg, cells, props, swat_init_, pressures);

                        const int np = props.numPhases();
                        for (int p = 0; p < np; ++p) {
                            copyFromRegion(pressures[p], cells, pp_[p]);
                            copyFromRegion(sat[p], cells, sat_[p]);
                        }
                        if (props.phaseUsage().phase_used[BlackoilPhases::Liquid]
                            && props.phaseUsage().phase_used[BlackoilPhases::Vapour]) {
                            const int oilpos = props.phaseUsage().phase_pos[BlackoilPhases::Liquid];
                            const int gaspos = props.phaseUsage().phase_pos[BlackoilPhases::Vapour];
                            const Vec rs_vals = computeRs(G, cells, pressures[oilpos], temp, *(rs_func_[r]), sat[gaspos]);
                            const Vec rv_vals = computeRs(G, cells, pressures[gaspos], temp, *(rv_func_[r]), sat[oilpos]);
                            copyFromRegion(rs_vals, cells, rs_);
                            copyFromRegion(rv_vals, cells, rv_);
                        }
                    }
                }

                template <class CellRangeType>
                void copyFromRegion(const Vec& source,
                                    const CellRangeType& cells,
                                    Vec& destination)
                {
                    auto s = source.begin();
                    auto c = cells.begin();
                    const auto e = cells.end();
                    for (; c != e; ++c, ++s) {
                        destination[*c] = *s;
                    }
                }

            };
        } // namespace DeckDependent
    } // namespace EQUIL
} // namespace Opm

#include <opm/core/simulator/initStateEquil_impl.hpp>

#endif // OPM_INITSTATEEQUIL_HEADER_INCLUDED
