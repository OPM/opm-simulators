/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FLOWDIAGNOSTICS_HEADER_INCLUDED
#define OPM_FLOWDIAGNOSTICS_HEADER_INCLUDED


#include <vector>
#include <utility>
#include <tuple>

struct Wells;

namespace Opm
{

    /// \brief Compute flow-capacity/storage-capacity based on time-of-flight.
    ///
    /// The F-Phi curve is an analogue to the fractional flow curve in a 1D
    /// displacement. It can be used to compute other interesting diagnostic
    /// quantities such as the Lorenz coefficient. For a technical description
    /// see Shavali et al. (SPE 146446), Shook and Mitchell (SPE 124625).
    ///
    /// \param[in]  pv    pore volumes of each cell
    /// \param[in]  ftof  forward (time from injector) time-of-flight values for each cell
    /// \param[in]  rtof  reverse (time to producer) time-of-flight values for each cell
    /// \return           a pair of vectors, the first containing F (flow capacity) the second
    ///                   containing Phi (storage capacity).
    std::pair<std::vector<double>, std::vector<double>>
    computeFandPhi(const std::vector<double>& pv,
                   const std::vector<double>& ftof,
                   const std::vector<double>& rtof);


    /// \brief Compute the Lorenz coefficient based on the F-Phi curve.
    ///
    /// The Lorenz coefficient is a measure of heterogeneity. It is equal
    /// to twice the area between the F-Phi curve and the F = Phi line.
    /// The coefficient can vary from zero to one. If the coefficient is
    /// zero (so the F-Phi curve is a straight line) we have perfect
    /// piston-like displacement while a coefficient of one indicates
    /// infinitely heterogenous displacement (essentially no sweep).
    ///
    /// Note: The coefficient is analogous to the Gini coefficient of
    /// economic theory, where the name Lorenz curve is applied to
    /// what we call the F-Phi curve.
    ///
    /// \param[in]  flowcap     flow capacity (F) as from computeFandPhi()
    /// \param[in]  storagecap  storage capacity (Phi) as from computeFandPhi()
    /// \return                 the Lorenz coefficient
    double computeLorenz(const std::vector<double>& flowcap,
                         const std::vector<double>& storagecap);


    /// \brief Compute sweep efficiency versus dimensionless time (PVI).
    ///
    /// The sweep efficiency is analogue to 1D displacement using the
    /// F-Phi curve as flux function.
    ///
    /// \param[in]  flowcap     flow capacity (F) as from computeFandPhi()
    /// \param[in]  storagecap  storage capacity (Phi) as from computeFandPhi()
    /// \return                 a pair of vectors, the first containing Ev (sweep efficiency)
    ///                         the second containing tD (dimensionless time).
    std::pair<std::vector<double>, std::vector<double>>
    computeSweep(const std::vector<double>& flowcap,
                 const std::vector<double>& storagecap);


    /// \brief Compute volumes associated with injector-producer pairs.
    ///
    /// \param[in]  wells       wells structure, containing NI injector wells and NP producer wells.
    /// \param[in]  porevol     pore volume of each grid cell
    /// \param[in]  ftracer     array of forward (injector) tracer values, NI per cell
    /// \param[in]  btracer     array of backward (producer) tracer values, NP per cell
    /// \return                 a vector of tuples, one tuple for each injector-producer pair,
    ///                         where the first and second elements are well indices for the
    ///                         injector and producer, and the third element is the pore volume
    ///                         associated with that pair.
    std::vector<std::tuple<int, int, double>>
    computeWellPairs(const Wells& wells,
                     const std::vector<double>& porevol,
                     const std::vector<double>& ftracer,
                     const std::vector<double>& btracer);

} // namespace Opm

#endif // OPM_FLOWDIAGNOSTICS_HEADER_INCLUDED
