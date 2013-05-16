/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#include <opm/autodiff/BlackoilPropsAd.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{

    // Making these typedef to make the code more readable.
    typedef BlackoilPropsAd::ADB ADB;
    typedef BlackoilPropsAd::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;

    /// Constructor wrapping an opm-core black oil interface.
    BlackoilPropsAd::BlackoilPropsAd(const BlackoilPropertiesInterface& props)
        : props_(props),
          pu_(props.phaseUsage())
    {
    }

    ////////////////////////////
    //      Rock interface    //
    ////////////////////////////

    /// \return   D, the number of spatial dimensions.
    int BlackoilPropsAd::numDimensions() const
    {
        return props_.numDimensions();
    }

    /// \return   N, the number of cells.
    int BlackoilPropsAd::numCells() const
    {
        return props_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropsAd::porosity() const
    {
        return props_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropsAd::permeability() const
    {
        return props_.permeability();
    }


    ////////////////////////////
    //      Fluid interface   //
    ////////////////////////////


    // ------ Density ------

    /// Densities of stock components at surface conditions.
    /// \return Array of 3 density values.
    const double* BlackoilPropsAd::surfaceDensity() const
    {
        return props_.surfaceDensity();
    }


    // ------ Viscosity ------

    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAd::muWat(const V& pw,
                             const Cells& cells) const
    {
        if (!pu_.phase_used[Water]) {
            THROW("Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        ASSERT(pw.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block mu(n, np);
        props_.viscosity(n, pw.data(), z.data(), cells.data(), mu.data(), 0);
        return mu.col(pu_.phase_pos[Water]);
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAd::muOil(const V& po,
                             const V& rs,
                             const Cells& cells) const
    {
        if (!pu_.phase_used[Oil]) {
            THROW("Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        ASSERT(po.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        if (pu_.phase_used[Gas]) {
            // Faking a z with the right ratio:
            //   rs = zg/zo
            z.col(pu_.phase_pos[Oil]) = V::Ones(n, 1);
            z.col(pu_.phase_pos[Gas]) = rs;
        }
        Block mu(n, np);
        props_.viscosity(n, po.data(), z.data(), cells.data(), mu.data(), 0);
        return mu.col(pu_.phase_pos[Oil]);
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    V BlackoilPropsAd::muGas(const V& pg,
                             const Cells& cells) const
    {
        if (!pu_.phase_used[Gas]) {
            THROW("Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        ASSERT(pg.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block mu(n, np);
        props_.viscosity(n, pg.data(), z.data(), cells.data(), mu.data(), 0);
        return mu.col(pu_.phase_pos[Gas]);
    }

    /// Water viscosity.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAd::muWat(const ADB& pw,
                               const Cells& cells) const
    {
        if (!pu_.phase_used[Water]) {
            THROW("Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        ASSERT(pw.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block mu(n, np);
        Block dmu(n, np);
        props_.viscosity(n, pw.value().data(), z.data(), cells.data(), mu.data(), dmu.data());
        ADB::M dmu_diag = spdiag(dmu.col(pu_.phase_pos[Water]));
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmu_diag * pw.derivative()[block];
        }
        return ADB::function(mu.col(pu_.phase_pos[Water]), jacs);
    }

    /// Oil viscosity.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAd::muOil(const ADB& po,
                               const ADB& rs,
                               const Cells& cells) const
    {
        if (!pu_.phase_used[Oil]) {
            THROW("Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        ASSERT(po.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        if (pu_.phase_used[Gas]) {
            // Faking a z with the right ratio:
            //   rs = zg/zo
            z.col(pu_.phase_pos[Oil]) = V::Ones(n, 1);
            z.col(pu_.phase_pos[Gas]) = rs.value();
        }
        Block mu(n, np);
        Block dmu(n, np);
        props_.viscosity(n, po.value().data(), z.data(), cells.data(), mu.data(), dmu.data());
        ADB::M dmu_diag = spdiag(dmu.col(pu_.phase_pos[Oil]));
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            // For now, we deliberately ignore the derivative with respect to rs,
            // since the BlackoilPropertiesInterface class does not evaluate it.
            // We would add to the next line: + dmu_drs_diag * rs.derivative()[block]
            jacs[block] = dmu_diag * po.derivative()[block];
        }
        return ADB::function(mu.col(pu_.phase_pos[Oil]), jacs);
    }

    /// Gas viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB BlackoilPropsAd::muGas(const ADB& pg,
                               const Cells& cells) const
    {
        if (!pu_.phase_used[Gas]) {
            THROW("Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        ASSERT(pg.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block mu(n, np);
        Block dmu(n, np);
        props_.viscosity(n, pg.value().data(), z.data(), cells.data(), mu.data(), dmu.data());
        ADB::M dmu_diag = spdiag(dmu.col(pu_.phase_pos[Gas]));
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = dmu_diag * pg.derivative()[block];
        }
        return ADB::function(mu.col(pu_.phase_pos[Gas]), jacs);
    }


    // ------ Formation volume factor (b) ------

    // These methods all call the matrix() method, after which the variable
    // (also) called 'matrix' contains, in each row, the A = RB^{-1} matrix for
    // a cell. For three-phase black oil:
    //  A = [  bw       0       0
    //          0       bo      0
    //          0      b0*rs   bw ]
    // Where b = B^{-1}.
    // Therefore, we extract the correct diagonal element, and are done.
    // When we need the derivatives (w.r.t. p, since we don't do w.r.t. rs),
    // we also get the following derivative matrix:
    //  A = [  dbw       0       0
    //          0       dbo      0
    //          0      db0*rs   dbw ]
    // Again, we just extract a diagonal element.

    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAd::bWat(const V& pw,
                            const Cells& cells) const
    {
        if (!pu_.phase_used[Water]) {
            THROW("Cannot call bWat(): water phase not present.");
        }
        const int n = cells.size();
        ASSERT(pw.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block matrix(n, np*np);
        props_.matrix(n, pw.data(), z.data(), cells.data(), matrix.data(), 0);
        const int wi = pu_.phase_pos[Water];
        return matrix.col(wi*np + wi);
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAd::bOil(const V& po,
                            const V& rs,
                            const Cells& cells) const
    {
        if (!pu_.phase_used[Oil]) {
            THROW("Cannot call bOil(): oil phase not present.");
        }
        const int n = cells.size();
        ASSERT(po.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        if (pu_.phase_used[Gas]) {
            // Faking a z with the right ratio:
            //   rs = zg/zo
            z.col(pu_.phase_pos[Oil]) = V::Ones(n, 1);
            z.col(pu_.phase_pos[Gas]) = rs;
        }
        Block matrix(n, np*np);
        props_.matrix(n, po.data(), z.data(), cells.data(), matrix.data(), 0);
        const int oi = pu_.phase_pos[Oil];
        return matrix.col(oi*np + oi);
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    V BlackoilPropsAd::bGas(const V& pg,
                            const Cells& cells) const
    {
        if (!pu_.phase_used[Gas]) {
            THROW("Cannot call bGas(): gas phase not present.");
        }
        const int n = cells.size();
        ASSERT(pg.size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block matrix(n, np*np);
        props_.matrix(n, pg.data(), z.data(), cells.data(), matrix.data(), 0);
        const int gi = pu_.phase_pos[Gas];
        return matrix.col(gi*np + gi);
    }

    /// Water formation volume factor.
    /// \param[in]  pw     Array of n water pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAd::bWat(const ADB& pw,
                              const Cells& cells) const
    {
        if (!pu_.phase_used[Water]) {
            THROW("Cannot call muWat(): water phase not present.");
        }
        const int n = cells.size();
        ASSERT(pw.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block matrix(n, np*np);
        Block dmatrix(n, np*np);
        props_.matrix(n, pw.value().data(), z.data(), cells.data(), matrix.data(), dmatrix.data());
        const int phase_ind = pu_.phase_pos[Water];
        const int column = phase_ind*np + phase_ind; // Index of our sought diagonal column.
        ADB::M db_diag = spdiag(dmatrix.col(column));
        const int num_blocks = pw.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = db_diag * pw.derivative()[block];
        }
        return ADB::function(matrix.col(column), jacs);
    }

    /// Oil formation volume factor.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  rs     Array of n gas solution factor values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAd::bOil(const ADB& po,
                              const ADB& rs,
                              const Cells& cells) const
    {
        if (!pu_.phase_used[Oil]) {
            THROW("Cannot call muOil(): oil phase not present.");
        }
        const int n = cells.size();
        ASSERT(po.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        if (pu_.phase_used[Gas]) {
            // Faking a z with the right ratio:
            //   rs = zg/zo
            z.col(pu_.phase_pos[Oil]) = V::Ones(n, 1);
            z.col(pu_.phase_pos[Gas]) = rs.value();
        }
        Block matrix(n, np*np);
        Block dmatrix(n, np*np);
        props_.matrix(n, po.value().data(), z.data(), cells.data(), matrix.data(), dmatrix.data());
        const int phase_ind = pu_.phase_pos[Oil];
        const int column = phase_ind*np + phase_ind; // Index of our sought diagonal column.
        ADB::M db_diag = spdiag(dmatrix.col(column));
        const int num_blocks = po.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            // For now, we deliberately ignore the derivative with respect to rs,
            // since the BlackoilPropertiesInterface class does not evaluate it.
            // We would add to the next line: + db_drs_diag * rs.derivative()[block]
            jacs[block] = db_diag * po.derivative()[block];
        }
        return ADB::function(matrix.col(column), jacs);
    }

    /// Gas formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB BlackoilPropsAd::bGas(const ADB& pg,
                              const Cells& cells) const
    {
        if (!pu_.phase_used[Gas]) {
            THROW("Cannot call muGas(): gas phase not present.");
        }
        const int n = cells.size();
        ASSERT(pg.value().size() == n);
        const int np = props_.numPhases();
        Block z = Block::Zero(n, np);
        Block matrix(n, np*np);
        Block dmatrix(n, np*np);
        props_.matrix(n, pg.value().data(), z.data(), cells.data(), matrix.data(), dmatrix.data());
        const int phase_ind = pu_.phase_pos[Gas];
        const int column = phase_ind*np + phase_ind; // Index of our sought diagonal column.
        ADB::M db_diag = spdiag(dmatrix.col(column));
        const int num_blocks = pg.numBlocks();
        std::vector<ADB::M> jacs(num_blocks);
        for (int block = 0; block < num_blocks; ++block) {
            jacs[block] = db_diag * pg.derivative()[block];
        }
        return ADB::function(matrix.col(column), jacs);
    }


#if 0
    // ------ Rs bubble point curve ------

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    V BlackoilPropsAd::rsMax(const V& po,
                             const Cells& cells) const
    {
    }

    /// Bubble point curve for Rs as function of oil pressure.
    /// \param[in]  po     Array of n oil pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n bubble point values for Rs.
    ADB BlackoilPropsAd::rsMax(const ADB& po,
                               const Cells& cells) const
    {
    }
#endif

    // ------ Relative permeability ------

    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  sg     Array of n gas saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 3 elements, each an array of n relperm values,
    ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
    std::vector<V> BlackoilPropsAd::relperm(const V& sw,
                                            const V& so,
                                            const V& sg,
                                            const Cells& cells) const
    {
        const int n = cells.size();
        const int np = props_.numPhases();
        Block s_all(n, np);
        if (pu_.phase_used[Water]) {
            ASSERT(sw.size() == n);
            s_all.col(pu_.phase_pos[Water]) = sw;
        }
        if (pu_.phase_used[Oil]) {
            ASSERT(so.size() == n);
            s_all.col(pu_.phase_pos[Oil]) = so;
        }
        if (pu_.phase_used[Gas]) {
            ASSERT(sg.size() == n);
            s_all.col(pu_.phase_pos[Gas]) = sg;
        }
        Block kr(n, np);
        props_.relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        std::vector<V> relperms;
        relperms.reserve(3);
        for (int phase = 0; phase < 3; ++phase) {
            if (pu_.phase_used[phase]) {
                relperms.emplace_back(kr.col(pu_.phase_pos[phase]));
            } else {
                relperms.emplace_back();
            }
        }
        return relperms;
    }

    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  sg     Array of n gas saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 3 elements, each an array of n relperm values,
    ///                    containing krw, kro, krg. Use PhaseIndex for indexing into the result.
    std::vector<ADB> BlackoilPropsAd::relperm(const ADB& sw,
                                              const ADB& so,
                                              const ADB& sg,
                                              const Cells& cells) const
    {
        const int n = cells.size();
        const int np = props_.numPhases();
        Block s_all(n, np);
        if (pu_.phase_used[Water]) {
            ASSERT(sw.value().size() == n);
            s_all.col(pu_.phase_pos[Water]) = sw.value();
        }
        if (pu_.phase_used[Oil]) {
            ASSERT(so.value().size() == n);
            s_all.col(pu_.phase_pos[Oil]) = so.value();
        } else {
            THROW("BlackoilPropsAd::relperm() assumes oil phase is active.");
        }
        if (pu_.phase_used[Gas]) {
            ASSERT(sg.value().size() == n);
            s_all.col(pu_.phase_pos[Gas]) = sg.value();
        }
        Block kr(n, np);
        Block dkr(n, np*np);
        props_.relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(3);
        typedef const ADB* ADBPtr;
        ADBPtr s[3] = { &sw, &so, &sg };
        for (int phase1 = 0; phase1 < 3; ++phase1) {
            if (pu_.phase_used[phase1]) {
                const int phase1_pos = pu_.phase_pos[phase1];
                std::vector<ADB::M> jacs(num_blocks);
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
                }
                for (int phase2 = 0; phase2 < 3; ++phase2) {
                    if (!pu_.phase_used[phase2]) {
                        continue;
                    }
                    const int phase2_pos = pu_.phase_pos[phase2];
                    // Assemble dkr1/ds2.
                    const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                    ADB::M dkr1_ds2_diag = spdiag(dkr.col(column));
                    for (int block = 0; block < num_blocks; ++block) {
                        jacs[block] += dkr1_ds2_diag * s[phase2]->derivative()[block];
                    }
                }
                relperms.emplace_back(ADB::function(kr.col(phase1_pos), jacs));
            } else {
                relperms.emplace_back(ADB::null());
            }
        }
        return relperms;
    }

} // namespace Opm

