/*
  Copyright 2015 IRIS
  
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

#ifndef SOLVENTPROPSADFROMDECK_HPP
#define SOLVENTPROPSADFROMDECK_HPP

#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>

#include <opm/core/props/pvt/PvtDead.hpp>
#include <opm/core/props/pvt/PvtInterface.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <cmath>
#include <vector>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{
class SolventPropsAdFromDeck
{
public:
    SolventPropsAdFromDeck(DeckConstPtr deck,
                           EclipseStateConstPtr eclipseState,
                           const int number_of_cells,
                           const int* global_cell);

    ////////////////////////////
    //      Fluid interface   //
    ////////////////////////////

    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef std::vector<int> Cells;
    typedef AutoDiffDenseBlock<double, 3> ADDB;

    /// Solvent formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADB bSolvent(const ADB& pg,
             const Cells& cells) const;

    /// Solvent viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADB muSolvent(const ADB& pg,
              const Cells& cells) const;

    /// Gas relPerm multipliers
    /// \param[in]  solventFraction Array of n solvent fraction values.
    /// \param[in]  cells           Array of n cell indices to be associated with the fraction values.
    /// \return                     Array of n gas relPerm multiplier values.
    ADB gasRelPermMultiplier(const ADB& solventFraction,
              const Cells& cells) const;

    /// Solvent relPerm multipliers
    /// \param[in]  solventFraction Array of n solvent fraction values.
    /// \param[in]  cells           Array of n cell indices to be associated with the fraction values.
    /// \return                     Array of n solvent relPerm multiplier values.
    ADB solventRelPermMultiplier(const ADB& solventFraction,
              const Cells& cells) const;


    // Dense block interfaces.

    /// Solvent formation volume factor.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n formation volume factor values.
    ADDB bSolvent(const ADDB& pg,
             const Cells& cells) const;

    /// Solvent viscosity.
    /// \param[in]  pg     Array of n gas pressure values.
    /// \param[in]  cells  Array of n cell indices to be associated with the pressure values.
    /// \return            Array of n viscosity values.
    ADDB muSolvent(const ADDB& pg,
              const Cells& cells) const;

    /// Gas relPerm multipliers
    /// \param[in]  solventFraction Array of n solvent fraction values.
    /// \param[in]  cells           Array of n cell indices to be associated with the fraction values.
    /// \return                     Array of n gas relPerm multiplier values.
    ADDB gasRelPermMultiplier(const ADDB& solventFraction,
              const Cells& cells) const;

    /// Solvent relPerm multipliers
    /// \param[in]  solventFraction Array of n solvent fraction values.
    /// \param[in]  cells           Array of n cell indices to be associated with the fraction values.
    /// \return                     Array of n solvent relPerm multiplier values.
    ADDB solventRelPermMultiplier(const ADDB& solventFraction,
              const Cells& cells) const;


    /// Solvent surface density
    /// \param[in]  cells           Array of n cell indices to be associated with the fraction values.
    /// \return                     Array of n solvent density values.
    V solventSurfaceDensity(const Cells& cells) const;

private:
    // The PVT region which is to be used for each cell
    std::vector<int> cellPvtRegionIdx_;
    std::vector<NonuniformTableLinear<double> > b_;
    std::vector<NonuniformTableLinear<double> > viscosity_;
    std::vector<NonuniformTableLinear<double> > inverseBmu_;
    std::vector<double> solvent_surface_densities_;
    std::vector<NonuniformTableLinear<double> > krg_;
    std::vector<NonuniformTableLinear<double> > krs_;
};

} // namespace OPM

#endif // SOLVENTPROPSADFROMDECK_HPP
