/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL.

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

#include <opm/autodiff/twophase/IncompPropsAdFromDeck.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{	
	/// Constructor wrapping an opm-core two-phase interface.
    IncompPropsAdFromDeck::IncompPropsAdFromDeck(const EclipseGridParser& deck,
                                                 const UnstructuredGrid& grid)
    {
        rock_.init(deck, grid);
        pvt_.init(deck);
        satprops_.init(deck, grid, 200);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "BlackoilPropsAdFromDeck::BlackoilPropsAdFromDeck() - "
                  "Inconsistent number of phases in pvt data (" << pvt_.numPhases()
                  << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
    }
    IncompPropsAdFromDeck::~IncompPropsAdFromDeck()
    {
    }

    ////////////////////////////
    //      Rock interface    //
    ////////////////////////////

    /// \return   D, the number of spatial dimensions.
    int IncompPropsAdFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int IncompPropsAdFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* IncompPropsAdFromDeck::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* IncompPropsAdFromDeck::permeability() const
    {
        return rock_.permeability();
    }

    ////////////////////////////
    //      Fluid interface   //
    ////////////////////////////
    
    /// \return   P, Number of active phases (also the number of components).
    int IncompPropsAdFromDeck::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return   Array of P viscosity values.
    const double* IncompPropsAdFromDeck::viscosity() const
    {
        return pvt_.viscosity();
    }

	///Densities of fluid phases at reservoir conditions.
    /// \return   Array of P density values.
    const double* IncompPropsAdFromDeck::density() const
    {
        return pvt_.reservoirDensities();
    }

	/// Densities of fluid phases at surface conditions.
    /// \return   Array of P density values.
    const double* IncompPropsAdFromDeck::surfaceDensity() const
    {
        return pvt_.surfaceDensities();
    }


    typedef IncompPropsAdFromDeck::ADB ADB;
    typedef IncompPropsAdFromDeck::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;


    // ------ Relative permeability ------
        
    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 2 elements, each an array of n relperm values,
    ///                    containing krw, kro.
    std::vector<V> 
    IncompPropsAdFromDeck::relperm(const V& sw,
                                   const V& so,
                                   const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw;
        s_all.col(1) = so;
        Block kr(n, np);
        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        std::vector<V> relperms;
        relperms.reserve(np);
        for (int phase = 0; phase < np; ++phase) {
            relperms.emplace_back(kr.col(phase));
        }
        return relperms;
    }



    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 2 elements, each an array of n relperm values,
    ///                    containing krw, kro.
    std::vector<ADB> 
    IncompPropsAdFromDeck::relperm(const ADB& sw,
                                   const ADB& so,
                                   const Cells& cells) const
    {
        const int n = cells.size();
        const int np = numPhases();
        Block s_all(n, np);
        assert(sw.size() == n && so.size() == n);
        s_all.col(0) = sw.value();
        s_all.col(1) = so.value();
        Block kr(n, np);
        Block dkr(n, np*np);
        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(np);
        typedef const ADB* ADBPtr;
        ADBPtr s[2] = { &sw, &so };
        for (int phase1 = 0; phase1 < np; ++phase1) {
            const int phase1_pos = phase1;
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
            }
            for (int phase2 = 0; phase2 < np; ++phase2) {
                const int phase2_pos = phase2;
                // Assemble dkr1/ds2.
                const int column = phase1_pos + np*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                ADB::M dkr1_ds2_diag = spdiag(dkr.col(column));
                for (int block = 0; block < num_blocks; ++block) {
                    jacs[block] += dkr1_ds2_diag * s[phase2]->derivative()[block];
                }
            }
                relperms.emplace_back(ADB::function(kr.col(phase1_pos), jacs));
        }
        return relperms;
    }

    /// Capillary pressure for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 2 elements, each an array of n capillary pressure values,
    ///                    containing the offsets for each p_o, p_w. The capillary pressure between
    ///                    two arbitrary phases alpha and beta is then given as p_alpha - p_beta.
    std::vector<ADB> 
	IncompPropsAdFromDeck::capPress(const ADB& sw,
                                    const ADB& so,
                                    const Cells& cells) const
    {
        const int numCells = cells.size();
        const int numActivePhases = numPhases();
        const int numBlocks = so.numBlocks();
        assert(sw.value().size() == numCells);
        assert(so.value().size() == numCells);
		Block s_all(numCells, numActivePhases);
        s_all.col(0) = sw.value();
        s_all.col(1) = so.value();

        Block pc(numCells, numActivePhases);
        Block dpc(numCells, numActivePhases*numActivePhases);
        satprops_.capPress(numCells, s_all.data(), cells.data(), pc.data(), dpc.data());

        std::vector<ADB> adbCapPressures;
        adbCapPressures.reserve(2);
        const ADB* s[2] = { &sw, &so};
        for (int phase1 = 0; phase1 < 2; ++phase1) {
            const int phase1_pos = phase1;
            std::vector<ADB::M> jacs(numBlocks);
            for (int block = 0; block < numBlocks; ++block) {
                jacs[block] = ADB::M(numCells, s[phase1]->derivative()[block].cols());
            }
            for (int phase2 = 0; phase2 < 2; ++phase2) {
                const int phase2_pos = phase2;
                    // Assemble dpc1/ds2.
                const int column = phase1_pos + numActivePhases*phase2_pos; // Recall: Fortran ordering from props_.relperm()
                ADB::M dpc1_ds2_diag = spdiag(dpc.col(column));
                for (int block = 0; block < numBlocks; ++block) {
                    jacs[block] += dpc1_ds2_diag * s[phase2]->derivative()[block];
                }
            }
            adbCapPressures.emplace_back(ADB::function(pc.col(phase1_pos), jacs));
        }
        return adbCapPressures;
    }


} //namespace Opm

