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

#include <opm/core/flowdiagnostics/FlowDiagnostics.hpp>
#include <opm/core/wells.h>

#include <opm/common/ErrorMacros.hpp>
#include <algorithm>
#include <numeric>

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
    std::pair<std::vector<double>, std::vector<double>> computeFandPhi(const std::vector<double>& pv,
                                                                       const std::vector<double>& ftof,
                                                                       const std::vector<double>& rtof)
    {
        if (pv.size() != ftof.size() || pv.size() != rtof.size()) {
            OPM_THROW(std::runtime_error, "computeFandPhi(): Input vectors must have same size.");
        }

        // Sort according to total travel time.
        const int n = pv.size();
        typedef std::pair<double, double> D2;
        std::vector<D2> time_and_pv(n);
        for (int ii = 0; ii < n; ++ii) {
            time_and_pv[ii].first = ftof[ii] + rtof[ii]; // Total travel time.
            time_and_pv[ii].second = pv[ii];
        }
        std::sort(time_and_pv.begin(), time_and_pv.end());

        // Compute Phi.
        std::vector<double> Phi(n + 1);
        Phi[0] = 0.0;
        for (int ii = 0; ii < n; ++ii) {
            Phi[ii+1] = time_and_pv[ii].second;
        }
        std::partial_sum(Phi.begin(), Phi.end(), Phi.begin());
        const double vt = Phi.back(); // Total pore volume.
        for (int ii = 1; ii < n+1; ++ii) { // Note limits of loop.
            Phi[ii] /= vt; // Normalize Phi.
        }

        // Compute F.
        std::vector<double> F(n + 1);
        F[0] = 0.0;
        for (int ii = 0; ii < n; ++ii) {
            F[ii+1] = time_and_pv[ii].second / time_and_pv[ii].first;
        }
        std::partial_sum(F.begin(), F.end(), F.begin());
        const double ft = F.back(); // Total flux.
        for (int ii = 1; ii < n+1; ++ii) { // Note limits of loop.
            F[ii] /= ft; // Normalize Phi.
        }

        return std::make_pair(F, Phi);
    }





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
                         const std::vector<double>& storagecap)
    {
        if (flowcap.size() != storagecap.size()) {
            OPM_THROW(std::runtime_error, "computeLorenz(): Input vectors must have same size.");
        }
        double integral = 0.0;
        // Trapezoid quadrature of the curve F(Phi).
        const int num_intervals = flowcap.size() - 1;
        for (int ii = 0; ii < num_intervals; ++ii) {
            const double len = storagecap[ii+1] - storagecap[ii];
            integral += (flowcap[ii] + flowcap[ii+1]) * len / 2.0;
        }
        return 2.0 * (integral - 0.5);
    }





    /// \brief Compute sweep efficiency versus dimensionless time (PVI).
    ///
    /// The sweep efficiency is analogue to 1D displacement using the
    /// F-Phi curve as flux function.
    ///
    /// \param[in]  flowcap     flow capacity (F) as from computeFandPhi()
    /// \param[in]  storagecap  storage capacity (Phi) as from computeFandPhi()
    /// \return                 a pair of vectors, the first containing Ev (sweep efficiency)
    ///                         the second containing tD (dimensionless time).
    std::pair<std::vector<double>, std::vector<double>> computeSweep(const std::vector<double>& flowcap,
                                                                     const std::vector<double>& storagecap)
    {
        if (flowcap.size() != storagecap.size()) {
            OPM_THROW(std::runtime_error, "computeSweep(): Input vectors must have same size.");
        }

        // Compute tD and Ev simultaneously,
        // skipping identical Phi data points.
        const int n = flowcap.size();
        std::vector<double> Ev;
        std::vector<double> tD;
        tD.reserve(n);
        Ev.reserve(n);
        tD.push_back(0.0);
        Ev.push_back(0.0);
        for (int ii = 1; ii < n; ++ii) { // Note loop limits.
            const double fd = flowcap[ii] - flowcap[ii-1];
            const double sd = storagecap[ii] - storagecap[ii-1];
            if (fd != 0.0) {
                tD.push_back(sd/fd);
                Ev.push_back(storagecap[ii] + (1.0 - flowcap[ii]) * tD.back());
            }
        }

        return std::make_pair(Ev, tD);
    }





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
    std::vector<std::tuple<int, int, double> >
    computeWellPairs(const Wells& wells,
                     const std::vector<double>& porevol,
                     const std::vector<double>& ftracer,
                     const std::vector<double>& btracer)
    {
        // Identify injectors and producers.
        std::vector<int> inj;
        std::vector<int> prod;
        const int nw = wells.number_of_wells;
        for (int w = 0; w < nw; ++w) {
            if (wells.type[w] == INJECTOR) {
                inj.push_back(w);
            } else {
                prod.push_back(w);
            }
        }

        // Check sizes of input arrays.
        const int nc = porevol.size();
        if (nc * inj.size() != ftracer.size()) {
            OPM_THROW(std::runtime_error, "computeWellPairs(): wrong size of input array ftracer.");
        }
        if (nc * prod.size() != btracer.size()) {
            OPM_THROW(std::runtime_error, "computeWellPairs(): wrong size of input array btracer.");
        }

        // Compute associated pore volumes.
        std::vector<std::tuple<int, int, double> > result;
        const int num_inj = inj.size();
        const int num_prod = prod.size();
        for (int inj_ix = 0; inj_ix < num_inj; ++inj_ix) {
            for (int prod_ix = 0; prod_ix < num_prod; ++prod_ix) {
                double assoc_porevol = 0.0;
                for (int c = 0; c < nc; ++c) {
                    assoc_porevol += porevol[c]
                        * ftracer[num_inj * c + inj_ix]
                        * btracer[num_prod * c + prod_ix];
                }
                result.push_back(std::make_tuple(inj[inj_ix], prod[prod_ix], assoc_porevol));
            }
        }
        return result;
    }



} // namespace Opm
