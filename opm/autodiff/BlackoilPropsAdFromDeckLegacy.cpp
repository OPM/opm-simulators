/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU.

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

#include <config.h>

#include <opm/autodiff/BlackoilPropsAdFromDeckLegacy.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>

#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/grid/utility/extractPvtTableIndex.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/common/ErrorMacros.hpp>

namespace Opm
{
    // Making these typedef to make the code more readable.
    typedef BlackoilPropsAdFromDeckLegacy::ADB ADB;
    typedef BlackoilPropsAdFromDeckLegacy::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeckLegacy::BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                                                 const Opm::EclipseState& eclState,
                                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                                 const UnstructuredGrid& grid,
                                                                 const bool init_rock)
      : BlackoilPropsAdFromDeck(deck,eclState,materialLawManager,grid,init_rock)
    {
    }

#ifdef HAVE_OPM_GRID
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeckLegacy::BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                                                 const Opm::EclipseState& eclState,
                                                                 const Dune::CpGrid& grid,
                                                                 const bool init_rock )
      : BlackoilPropsAdFromDeck(deck,eclState,grid,init_rock)
    {
    }
#endif

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeckLegacy::BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                                                 const Opm::EclipseState& eclState,
                                                                 const UnstructuredGrid& grid,
                                                                 const bool init_rock)
      : BlackoilPropsAdFromDeck(deck,eclState,grid,init_rock)
    {
    }

#ifdef HAVE_OPM_GRID
    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAdFromDeckLegacy::BlackoilPropsAdFromDeckLegacy(const Opm::Deck& deck,
                                                                 const Opm::EclipseState& eclState,
                                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                                 const Dune::CpGrid& grid,
                                                                 const bool init_rock )
      : BlackoilPropsAdFromDeck(deck,eclState,materialLawManager,grid,init_rock)
    {
    }
#endif

/// Constructor for properties on a subgrid
    BlackoilPropsAdFromDeckLegacy::BlackoilPropsAdFromDeckLegacy(const BlackoilPropsAdFromDeckLegacy& props,
                                                             std::shared_ptr<MaterialLawManager> materialLawManager,
                                                             const int number_of_cells)
    : BlackoilPropsAdFromDeck(props,materialLawManager,number_of_cells)
    {
    }

    // ------ Density ------

    /// Densities of stock components at surface conditions.
    /// \param[in] phaseIdx
    /// \param[in] cells  Array of n cell indices to be associated with the pressure values.
    /// \return Array of n density values for phase given by phaseIdx.
    V BlackoilPropsAdFromDeckLegacy::surfaceDensity(const int phaseIdx, const Cells& cells) const
    {
        assert( !(phaseIdx > numPhases()));
        const int n = cells.size();
        V rhos = V::Zero(n);
        for (int cellIdx = 0; cellIdx < n; ++cellIdx) {
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            rhos[cellIdx] = FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
        }
        return rhos;
    }


    // ------ Viscosity ------


    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeckLegacy::muWat(const ADB& pw,
                                             const ADB& T,
                                             const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call muWat(): water phase not active.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V mu(n);
        V dmudp(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/1> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;

        pEval.setDerivative(0, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(pw.value()[i]);
            TEval.setValue(T.value()[i]);

            const Eval& muEval = FluidSystem::waterPvt().viscosity(pvtRegionIdx, TEval, pEval);

            mu[i] = muEval.value();
            dmudp[i] = muEval.derivative(0);
        }

        if (pw.derivative().empty()) {
            return ADB::constant(std::move(mu));
        } else {
            ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
            const int num_blocks = pw.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                fastSparseProduct(dmudp_diag, pw.derivative()[block], jacs[block]);
            }
            return ADB::function(std::move(mu), std::move(jacs));
        }
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeckLegacy::muOil(const ADB& po,
                                             const ADB& T,
                                             const ADB& rs,
                                             const std::vector<PhasePresence>& cond,
                                             const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call muOil(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/2> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;
        Eval RsEval = 0.0;

        pEval.setDerivative(0, 1.0);
        RsEval.setDerivative(1, 1.0);

        Eval muEval;
        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(po.value()[i]);
            TEval.setValue(T.value()[i]);

            if (cond[i].hasFreeGas()) {
                muEval = FluidSystem::oilPvt().saturatedViscosity(pvtRegionIdx, TEval, pEval);
            }
            else {
                if (phase_usage_.phase_used[Gas]) {
                    RsEval.setValue(rs.value()[i]);
                }
                muEval = FluidSystem::oilPvt().viscosity(pvtRegionIdx, TEval, pEval, RsEval);
            }

            mu[i] = muEval.value();
            dmudp[i] = muEval.derivative(0);
            dmudr[i] = muEval.derivative(1);
        }

        ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
        ADB::M dmudr_diag(dmudr.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dmudp_diag, po.derivative()[block], jacs[block]);
            if (phase_usage_.phase_used[Gas]) {
                ADB::M temp;
                fastSparseProduct(dmudr_diag, rs.derivative()[block], temp);
                jacs[block] += temp;
            }
        }
        return ADB::function(std::move(mu), std::move(jacs));
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAdFromDeckLegacy::muGas(const ADB& pg,
                                             const ADB& T,
                                             const ADB& rv,
                                             const std::vector<PhasePresence>& cond,
                                             const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call muGas(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.value().size() == n);
        V mu(n);
        V dmudp(n);
        V dmudr(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/2> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;
        Eval RvEval = 0.0;
        Eval muEval;

        pEval.setDerivative(0, 1.0);
        RvEval.setDerivative(1, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(pg.value()[i]);
            TEval.setValue(T.value()[i]);

            if (cond[i].hasFreeOil()) {
                muEval = FluidSystem::gasPvt().saturatedViscosity(pvtRegionIdx, TEval, pEval);
            }
            else {
                RvEval.setValue(rv.value()[i]);
                muEval = FluidSystem::gasPvt().viscosity(pvtRegionIdx, TEval, pEval, RvEval);
            }

            mu[i] = muEval.value();
            dmudp[i] = muEval.derivative(0);
            dmudr[i] = muEval.derivative(1);
        }

        ADB::M dmudp_diag(dmudp.matrix().asDiagonal());
        ADB::M dmudr_diag(dmudr.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dmudp_diag, pg.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dmudr_diag, rv.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(mu), std::move(jacs));
    }


    // ------ Formation volume factor (b) ------


    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeckLegacy::bWat(const ADB& pw,
                                            const ADB& T,
                                            const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Water]) {
            OPM_THROW(std::runtime_error, "Cannot call bWat(): water phase not active.");
        }
        const int n = cells.size();
        assert(pw.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/1> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;

        pEval.setDerivative(0, 1.0);
        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(pw.value()[i]);
            TEval.setValue(T.value()[i]);

            const Eval& bEval = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, TEval, pEval);

            b[i] = bEval.value();
            dbdp[i] = bEval.derivative(0);
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, pw.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(b), std::move(jacs));
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cond   Array of n taxonomies classifying fluid condition.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeckLegacy::bOil(const ADB& po,
                                            const ADB& T,
                                            const ADB& rs,
                                            const std::vector<PhasePresence>& cond,
                                            const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call bOil(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/2> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;
        Eval RsEval = 0.0;
        Eval bEval;

        pEval.setDerivative(0, 1.0);
        RsEval.setDerivative(1, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(po.value()[i]);
            TEval.setValue(T.value()[i]);

            //RS/RV only makes sense when gas phase is active
            if (cond[i].hasFreeGas()) {
                bEval = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx, TEval, pEval);
            }
            else {
                if (rs.size() == 0) {
                    RsEval.setValue(0.0);
                }
                else {
                    RsEval.setValue(rs.value()[i]);
                }
                bEval = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, TEval, pEval, RsEval);
            }

            b[i] = bEval.value();
            dbdp[i] = bEval.derivative(0);
            dbdr[i] = bEval.derivative(1);
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        ADB::M dbdr_diag(dbdr.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);

        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, po.derivative()[block], jacs[block]);
            if (phase_usage_.phase_used[Gas]) {
                ADB::M temp;
                fastSparseProduct(dbdr_diag, rs.derivative()[block], temp);
                jacs[block] += temp;
            }
        }
        return ADB::function(std::move(b), std::move(jacs));
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  rv     Array of n vapor oil/gas ratio
    /// \param[in]  cond   Array of n objects, each specifying which phases are present with non-zero saturation in a cell.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAdFromDeckLegacy::bGas(const ADB& pg,
                                            const ADB& T,
                                            const ADB& rv,
                                            const std::vector<PhasePresence>& cond,
                                            const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call bGas(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.size() == n);

        V b(n);
        V dbdp(n);
        V dbdr(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/2> Eval;

        Eval pEval = 0.0;
        Eval TEval = 0.0;
        Eval RvEval = 0.0;
        Eval bEval;

        pEval.setDerivative(0, 1.0);
        RvEval.setDerivative(1, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(pg.value()[i]);
            TEval.setValue(T.value()[i]);

            if (cond[i].hasFreeOil()) {
                bEval = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvtRegionIdx, TEval, pEval);
            }
            else {
                RvEval.setValue(rv.value()[i]);
                bEval = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, TEval, pEval, RvEval);
            }

            b[i] = bEval.value();
            dbdp[i] = bEval.derivative(0);
            dbdr[i] = bEval.derivative(1);
        }

        ADB::M dbdp_diag(dbdp.matrix().asDiagonal());
        ADB::M dbdr_diag(dbdr.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(dbdp_diag, pg.derivative()[block], jacs[block]);
            ADB::M temp;
            fastSparseProduct(dbdr_diag, rv.derivative()[block], temp);
            jacs[block] += temp;
        }
        return ADB::function(std::move(b), std::move(jacs));
    }



    // ------ Rs bubble point curve ------

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAdFromDeckLegacy::rsSat(const ADB& po,
                                             const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Oil]) {
            OPM_THROW(std::runtime_error, "Cannot call rsSat(): oil phase not active.");
        }
        const int n = cells.size();
        assert(po.size() == n);
        V rbub(n);
        V drbubdp(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/1> Eval;

        Eval pEval = 0.0;
        Eval TEval = 293.15; // temperature is not supported by this API!

        pEval.setDerivative(0, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(po.value()[i]);

            const Eval& RsEval = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx, TEval, pEval);

            rbub[i] = RsEval.value();
            drbubdp[i] = RsEval.derivative(0);
        }

        ADB::M drbubdp_diag(drbubdp.matrix().asDiagonal());
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(drbubdp_diag, po.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(rbub), std::move(jacs));
    }

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAdFromDeckLegacy::rsSat(const ADB& po,
                                             const ADB& so,
                                             const Cells& cells) const
    {
        ADB rs = rsSat(po, cells);
        applyVap(rs, so, cells, vap2_);
        return rs;
    }

    // ------ Rv condensation curve ------

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n condensation point values for Rv.
    ADB BlackoilPropsAdFromDeckLegacy::rvSat(const ADB& pg,
                                             const Cells& cells) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call rvSat(): gas phase not active.");
        }
        const int n = cells.size();
        assert(pg.size() == n);
        V rv(n);
        V drvdp(n);

        typedef Opm::DenseAd::Evaluation<double, /*size=*/1> Eval;

        Eval pEval = 0.0;
        Eval TEval = 293.15; // temperature is not supported by this API!

        pEval.setDerivative(0, 1.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            pEval.setValue(pg.value()[i]);

            const Eval& RvEval = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx, TEval, pEval);

            rv[i] = RvEval.value();
            drvdp[i] = RvEval.derivative(0);
        }

        ADB::M drvdp_diag(drvdp.matrix().asDiagonal());
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(drvdp_diag, pg.derivative()[block], jacs[block]);
        }
        return ADB::function(std::move(rv), std::move(jacs));
    }

    /// Condensation curve for Rv as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n condensation point values for Rv.
    ADB BlackoilPropsAdFromDeckLegacy::rvSat(const ADB& po,
                                             const ADB& so,
                                             const Cells& cells) const
    {
        ADB rv = rvSat(po, cells);
        applyVap(rv, so, cells, vap1_);
        return rv;
    }

    // ------ Relative permeability ------

    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  sg     Array of n gas saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 3 elements, each an array of n relperm values,
    ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
    std::vector<ADB> BlackoilPropsAdFromDeckLegacy::relperm(const ADB& sw,
                                                            const ADB& so,
                                                            const ADB& sg,
                                                            const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Water]) = sw.value();
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Oil]) = so.value();
        } else {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::relperm() assumes oil phase is active.");
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.value().size() == n);
            s_all.col(phase_usage_.phase_pos[Gas]) = sg.value();
        }
        Block kr(n, np);
        Block dkr(n, np*np);
        satprops_->relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(3);
        typedef const ADB* ADBPtr;
        ADBPtr s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (phase_usage_.phase_used[phase1]) {
                const int phase1_pos = phase_usage_.phase_pos[phase1];
                std::vector<ADB::M> jacs(num_blocks);
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!phase_usage_.phase_used[phase2]) {
                        continue;
                    }
                    const int phase2_pos = phase_usage_.phase_pos[phase2];
                    // Assemble dkr1/ds2.
                    const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()

                    ADB::M dkr1_ds2_diag(dkr.col(column).matrix().asDiagonal());
                    for (int block = 0; block < num_blocks; ++block) {
                        ADB::M temp;
                        fastSparseProduct(dkr1_ds2_diag, s[phase2]->derivative()[block], temp);
                        jacs[block] += temp;
                    }
                }
                ADB::V val = kr.col(phase1_pos);
                relperms.emplace_back(ADB::function(std::move(val), std::move(jacs)));
            } else {
                relperms.emplace_back(ADB::null());
            }
        }
        return relperms;
    }

    std::vector<ADB> BlackoilPropsAdFromDeckLegacy::capPress(const ADB& sw,
                                                             const ADB& so,
                                                             const ADB& sg,
                                                             const Cells& cells) const
    {
        const int nCells = cells.size();
        const int nActivePhases = numPhases();
        const int nBlocks = so.numBlocks();

        Block activeSat(nCells, nActivePhases);
        if (phase_usage_.phase_used[Water]) {
            assert(sw.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Water]) = sw.value();
        }
        if (phase_usage_.phase_used[Oil]) {
            assert(so.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Oil]) = so.value();
        } else {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::relperm() assumes oil phase is active.");
        }
        if (phase_usage_.phase_used[Gas]) {
            assert(sg.value().size() == nCells);
            activeSat.col(phase_usage_.phase_pos[Gas]) = sg.value();
        }

        Block pc(nCells, nActivePhases);
        Block dpc(nCells, nActivePhases*nActivePhases);
        satprops_->capPress(nCells, activeSat.data(), cells.data(), pc.data(), dpc.data());

        std::vector<ADB> adbCapPressures;
        adbCapPressures.reserve(3);
        const ADB* s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (phase_usage_.phase_used[phase1]) {
                const int phase1_pos = phase_usage_.phase_pos[phase1];
                std::vector<ADB::M> jacs(nBlocks);
                for (int block = 0; block < nBlocks; ++block) {
                    jacs[block] = ADB::M(nCells, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!phase_usage_.phase_used[phase2])
                        continue;
                    const int phase2_pos = phase_usage_.phase_pos[phase2];
                    // Assemble dpc1/ds2.
                    const int column = phase1_pos + nActivePhases*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                    ADB::M dpc1_ds2_diag(dpc.col(column).matrix().asDiagonal());
                    for (int block = 0; block < nBlocks; ++block) {
                        ADB::M temp;
                        fastSparseProduct(dpc1_ds2_diag, s[phase2]->derivative()[block], temp);
                        jacs[block] += temp;
                    }
                }
                ADB::V val = pc.col(phase1_pos);
                adbCapPressures.emplace_back(ADB::function(std::move(val), std::move(jacs)));
            } else {
                adbCapPressures.emplace_back(ADB::null());
            }
        }
        return adbCapPressures;
    }

    /// Bubble point pressures
    std::vector<double>
    BlackoilPropsAdFromDeckLegacy::bubblePointPressure(const Cells& cells,
                                                       const V& T,
                                                       const V& rs) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call bubblePointPressure(): gas phase not active.");
        }
        const int n = cells.size();
        std::vector<double> Pb(n, 0.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            try {
                Pb[i] = FluidSystem::oilPvt().saturationPressure(pvtRegionIdx, T[i], rs[i]);
            }
            catch (const NumericalIssue&) {
                // Ignore
            }
        }

        return Pb;
    }

    /// Dew point pressures
    std::vector<double>
    BlackoilPropsAdFromDeckLegacy::dewPointPressure(const Cells& cells,
                                                    const V& T,
                                                    const V& rv) const
    {
        if (!phase_usage_.phase_used[Gas]) {
            OPM_THROW(std::runtime_error, "Cannot call dewPointPressure(): gas phase not active.");
        }
        const int n = cells.size();
        std::vector<double> Pd(n, 0.0);

        for (int i = 0; i < n; ++i) {
            unsigned pvtRegionIdx = cellPvtRegionIdx_[cells[i]];
            try {
                Pd[i] = FluidSystem::gasPvt().saturationPressure(pvtRegionIdx, T[i], rv[i]);
            }
            catch (const NumericalIssue&) {
                // Ignore
            }
        }

        return Pd;
    }

    /// Apply correction to rs/rv according to kw VAPPARS
    /// \param[in/out] r     Array of n rs/rv values.
    /// \param[in]     so    Array of n oil saturation values.
    /// \param[in]     cells Array of n cell indices to be associated with the r and so values.
    /// \param[in]     vap   Correction parameter.
    void BlackoilPropsAdFromDeckLegacy::applyVap(V& r,
                                                 const V& so,
                                                 const std::vector<int>& cells,
                                                 const double vap) const
    {
        if (vap > 0.0) {
            const int n = cells.size();
            V factor = V::Ones(n, 1);
            const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
            for (int i=0; i<n; ++i) {
                if (satOilMax_[cells[i]] > vap_satmax_guard_ && so[i] < satOilMax_[cells[i]]) {
                    // guard against too small saturation values.
                    const double so_i= std::max(so[i],eps_sqrt);
                    factor[i] = std::pow(so_i/satOilMax_[cells[i]], vap);
                }
            }
            r = factor*r;
        }
    }

    /// Apply correction to rs/rv according to kw VAPPARS
    /// \param[in/out] r     Array of n rs/rv values.
    /// \param[in]     so    Array of n oil saturation values.
    /// \param[in]     cells Array of n cell indices to be associated with the r and so values.
    /// \param[in]     vap   Correction parameter.
    void BlackoilPropsAdFromDeckLegacy::applyVap(ADB& r,
                                                 const ADB& so,
                                                 const std::vector<int>& cells,
                                                 const double vap) const
    {
        if (vap > 0.0) {
            const int n = cells.size();
            V factor = V::Ones(n, 1);
            const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
            V dfactor_dso = V::Zero(n, 1);
            for (int i=0; i<n; ++i) {
                if (satOilMax_[cells[i]] > vap_satmax_guard_ && so.value()[i] < satOilMax_[cells[i]]) {
                    // guard against too small saturation values.
                    const double so_i= std::max(so.value()[i],eps_sqrt);
                    factor[i] = std::pow(so_i/satOilMax_[cells[i]], vap);
                    dfactor_dso[i] = vap*std::pow(so_i/satOilMax_[cells[i]], vap-1.0)/satOilMax_[cells[i]];
                }
            }
            ADB::M dfactor_dso_diag(dfactor_dso.matrix().asDiagonal());
            const int num_blocks = so.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = dfactor_dso_diag * so.derivative()[block];
            }
            r = ADB::function(std::move(factor), std::move(jacs))*r;
        }
    }

    /// Obtain the scaled critical oil in gas saturation values.
    /// \param[in]  cells  Array of cell indices.
    /// \return Array of critical oil in gas saturaion values.
    V BlackoilPropsAdFromDeckLegacy::scaledCriticalOilinGasSaturations(const Cells& cells) const {
        assert(phase_usage_.phase_used[Gas]);
        assert(phase_usage_.phase_used[Oil]);

        const int n = cells.size();
        V sogcr = V::Zero(n);
        const MaterialLawManager& materialLawManager = satprops_->materialLawManager();
        for (int i = 0; i < n; ++i) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cells[i]);
                sogcr[i] = scaledDrainageInfo.Sogcr;
        }
        return sogcr;
    }

    /// Obtain the scaled critical gas saturation values.
    /// \param[in]  cells  Array of cell indices.
    /// \return Array of scaled critical gas saturaion values.
    V BlackoilPropsAdFromDeckLegacy::scaledCriticalGasSaturations(const Cells& cells) const {
        assert(phase_usage_.phase_used[Gas]);

        const int n = cells.size();
        V sgcr = V::Zero(n);
        const MaterialLawManager& materialLawManager = satprops_->materialLawManager();
        for (int i = 0; i < n; ++i) {
            const auto& scaledDrainageInfo =
                materialLawManager.oilWaterScaledEpsInfoDrainage(cells[i]);
                sgcr[i] = scaledDrainageInfo.Sgcr;
        }
        return sgcr;
    }


} // namespace Opm
