/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#include <opm/autodiff/IncompPropsAdBasic.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>

namespace Opm
{
	/// Constructor.
    IncompPropsAdBasic::IncompPropsAdBasic(const parameter::ParameterGroup& param,
                                           const int dim,
                                           const int num_cells)
    {
        double poro = param.getDefault("porosity", 1.0);
        using namespace Opm::unit;
        using namespace Opm::prefix;
        double perm = param.getDefault("permeability", 100) * milli * darcy;
        rock_.init(dim, num_cells, poro, perm);
        pvt_.init(param);
        satprops_.init(param);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "IncompPropsAdBasic::IncompPropsAdBasic() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
        viscosity_.resize(pvt_.numPhases());
        pvt_.mu(1, 0, 0, &viscosity_[0]);
    }

	/// Constructor.
    IncompPropsAdBasic::IncompPropsAdBasic(const int num_phases,
                                                 const  SaturationPropsBasic::RelPermFunc& relpermfunc,
                                                 const std::vector<double>&  rho,
                                                 const std::vector<double>& mu,
                                                 const double por, //porosity
                                                 const double perm,
                                                 const int dim,
                                                 const int num_cells)
    {
        rock_.init(dim, num_cells, por, perm);
        pvt_.init(num_phases, rho, mu);
        satprops_.init(num_phases, relpermfunc);
        if (pvt_.numPhases() != satprops_.numPhases()) {
            OPM_THROW(std::runtime_error, "IncompPropsAdBasic::IncompPropsAdBasic() - Inconsistent number of phases in pvt data ("
                  << pvt_.numPhases() << ") and saturation-dependent function data (" << satprops_.numPhases() << ").");
        }
        viscosity_.resize(pvt_.numPhases());
        pvt_.mu(1, 0, 0, &viscosity_[0]);
    }

	/// Destructor.
    IncompPropsAdBasic::~IncompPropsAdBasic()
    {
    }

    ////////////////////////////
    //      Rock interface    //
    ////////////////////////////
    
    /// \return   D, the number of spatial dimensions.
    int IncompPropsAdBasic::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int IncompPropsAdBasic::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* IncompPropsAdBasic::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* IncompPropsAdBasic::permeability() const
    {
        return rock_.permeability();
    }


    ////////////////////////////
    //      Fluid interface   //
    ////////////////////////////

    /// \return   P, the number of phases (also the number of components).
    int IncompPropsAdBasic::numPhases() const
    {
        return pvt_.numPhases();
    }

    /// \return Array of P viscosity values.
    const double* IncompPropsAdBasic::viscosity() const
    {
        return &viscosity_[0];
    }

	/// Densities of fluid phases at reservoir conditions.
    /// \return   Array of P density values.
    const double* IncompPropsAdBasic::density() const
    {
        return pvt_.surfaceDensities();
    }

	/// Densities of fluid phases at surface conditions.
    /// \return   Array of P density values.
    const double* IncompPropsAdBasic::surfaceDensity() const
    {
        return pvt_.surfaceDensities();
    }

    typedef IncompPropsAdBasic::ADB ADB;
    typedef IncompPropsAdBasic::V V;
    typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Block;
    typedef std::vector<int> Cells;

    // ------ Relative permeability ------
        
    /// Relative permeabilities for all phases.
    /// \param[in]  sw     Array of n water saturation values.
    /// \param[in]  so     Array of n oil saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the saturation values.
    /// \return            An std::vector with 2 elements, each an array of n relperm values,
    ///                    containing krw, kro.
    std::vector<V> 
    IncompPropsAdBasic::relperm(const V& sw,
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
//        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), 0);
        satprops_.relperm(n, s_all.data(), kr.data(), 0);

        std::vector<V> relperms;
        relperms.reserve(2);
        for (int phase = 0; phase < 2; ++phase) {
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
    IncompPropsAdBasic::relperm(const ADB& sw,
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
//        satprops_.relperm(n, s_all.data(), cells.data(), kr.data(), dkr.data());
        satprops_.relperm(n, s_all.data(), kr.data(), dkr.data());
        const int num_blocks = so.numBlocks();
        std::vector<ADB> relperms;
        relperms.reserve(2);
        typedef const ADB* ADBPtr;
        ADBPtr s[2] = { &sw, &so };
        for (int phase1 = 0; phase1 < 2; ++phase1) {
            const int phase1_pos = phase1;
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = ADB::M(n, s[phase1]->derivative()[block].cols());
            }
            for (int phase2 = 0; phase2 < 2; ++phase2) {
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

} //namespace Opm
