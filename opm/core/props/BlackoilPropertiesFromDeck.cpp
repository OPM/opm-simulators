/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include "config.h"
#include <opm/core/props/BlackoilPropertiesFromDeck.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/compressedToCartesian.hpp>
#include <opm/core/utility/extractPvtTableIndex.hpp>
#include <vector>
#include <numeric>

namespace Opm
{
    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           const UnstructuredGrid& grid,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx
            = compressedToCartesian(grid.number_of_cells, grid.global_cell);

        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck, eclState, materialLawManager, grid.number_of_cells, grid.global_cell, grid.cartdims,
             init_rock);
    }

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           const UnstructuredGrid& grid,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx
            = compressedToCartesian(grid.number_of_cells, grid.global_cell);

        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck, eclState, materialLawManager, grid.number_of_cells, grid.global_cell, grid.cartdims, param, init_rock);
    }

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx
            = compressedToCartesian(number_of_cells, global_cell);

        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck, eclState, materialLawManager, number_of_cells, global_cell, cart_dims,
             init_rock);
    }

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        std::vector<int> compressedToCartesianIdx
            = compressedToCartesian(number_of_cells, global_cell);

        auto materialLawManager = std::make_shared<MaterialLawManager>();
        materialLawManager->initFromDeck(deck, eclState, compressedToCartesianIdx);

        init(deck,
             eclState,
             materialLawManager,
             number_of_cells,
             global_cell,
             cart_dims,
             param,
             init_rock);
    }

    BlackoilPropertiesFromDeck::BlackoilPropertiesFromDeck(Opm::DeckConstPtr deck,
                                                           Opm::EclipseStateConstPtr eclState,
                                                           std::shared_ptr<MaterialLawManager> materialLawManager,
                                                           int number_of_cells,
                                                           const int* global_cell,
                                                           const int* cart_dims,
                                                           const parameter::ParameterGroup& param,
                                                           bool init_rock)
    {
        init(deck,
             eclState,
             materialLawManager,
             number_of_cells,
             global_cell,
             cart_dims,
             param,
             init_rock);
    }

    inline void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr deck,
                                                 Opm::EclipseStateConstPtr eclState,
                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                 int number_of_cells,
                                                 const int* global_cell,
                                                 const int* cart_dims,
                                                 bool init_rock)
    {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, eclState, number_of_cells, global_cell);

        if (init_rock){
           rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }
        phaseUsage_ = phaseUsageFromDeck(deck);
        initSurfaceDensities_(deck);
        oilPvt_.initFromDeck(deck, eclState);
        gasPvt_.initFromDeck(deck, eclState);
        waterPvt_.initFromDeck(deck, eclState);
        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        ptr->init(phaseUsageFromDeck(deck), materialLawManager);
        satprops_.reset(ptr);
    }

    inline void BlackoilPropertiesFromDeck::init(Opm::DeckConstPtr deck,
                                                 Opm::EclipseStateConstPtr eclState,
                                                 std::shared_ptr<MaterialLawManager> materialLawManager,
                                                 int number_of_cells,
                                                 const int* global_cell,
                                                 const int* cart_dims,
                                                 const parameter::ParameterGroup& param,
                                                 bool init_rock)
    {
        // retrieve the cell specific PVT table index from the deck
        // and using the grid...
        extractPvtTableIndex(cellPvtRegionIdx_, eclState, number_of_cells, global_cell);

        if(init_rock){
            rock_.init(eclState, number_of_cells, global_cell, cart_dims);
        }

        phaseUsage_ = phaseUsageFromDeck(deck);
        initSurfaceDensities_(deck);
        oilPvt_.initFromDeck(deck, eclState);
        gasPvt_.initFromDeck(deck, eclState);
        waterPvt_.initFromDeck(deck, eclState);

        // Unfortunate lack of pointer smartness here...
        std::string threephase_model = param.getDefault<std::string>("threephase_model", "gwseg");
        if (deck->hasKeyword("ENDSCALE") && threephase_model != "gwseg") {
            OPM_THROW(std::runtime_error, "Sorry, end point scaling currently available for the 'gwseg' model only.");
        }

        SaturationPropsFromDeck* ptr
            = new SaturationPropsFromDeck();
        ptr->init(phaseUsageFromDeck(deck), materialLawManager);
        satprops_.reset(ptr);
    }

    BlackoilPropertiesFromDeck::~BlackoilPropertiesFromDeck()
    {
    }

    /// \return   D, the number of spatial dimensions.
    int BlackoilPropertiesFromDeck::numDimensions() const
    {
        return rock_.numDimensions();
    }

    /// \return   N, the number of cells.
    int BlackoilPropertiesFromDeck::numCells() const
    {
        return rock_.numCells();
    }

    /// \return   Array of N porosity values.
    const double* BlackoilPropertiesFromDeck::porosity() const
    {
        return rock_.porosity();
    }

    /// \return   Array of ND^2 permeability values.
    ///           The D^2 permeability values for a cell are organized as a matrix,
    ///           which is symmetric (so ordering does not matter).
    const double* BlackoilPropertiesFromDeck::permeability() const
    {
        return rock_.permeability();
    }


    // ---- Fluid interface ----

    /// \return   P, the number of phases (also the number of components).
    int BlackoilPropertiesFromDeck::numPhases() const
    {
        return phaseUsage_.num_phases;
    }

    /// \return   Object describing the active phases.
    PhaseUsage BlackoilPropertiesFromDeck::phaseUsage() const
    {
        return phaseUsage_;
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  p      Array of n pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  z      Array of nP surface volume values.
    /// \param[in]  cells  Array of n cell indices to be associated with the p and z values.
    /// \param[out] mu     Array of nP viscosity values, array must be valid before calling.
    /// \param[out] dmudp  If non-null: array of nP viscosity derivative values,
    ///                    array must be valid before calling.
    void BlackoilPropertiesFromDeck::viscosity(const int n,
                                               const double* p,
                                               const double* T,
                                               const double* z,
                                               const int* cells,
                                               double* mu,
                                               double* dmudp) const
    {
        const auto& pu = phaseUsage();

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RsLad = 0.0;
        LadEval RvLad = 0.0;
        LadEval muLad = 0.0;

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++ i) {
            int cellIdx = cells[i];
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            pLad.value = p[i];
            TLad.value = T[i];

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                muLad = waterPvt_.viscosity(pvtRegionIdx, TLad, pLad);
                int offset = pu.num_phases*cellIdx + pu.phase_pos[BlackoilPhases::Aqua];
                mu[offset] = muLad.value;
                dmudp[offset] = muLad.derivatives[0];
            }

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                muLad = oilPvt_.viscosity(pvtRegionIdx, TLad, pLad, RsLad);
                int offset = pu.num_phases*cellIdx + pu.phase_pos[BlackoilPhases::Liquid];
                mu[offset] = muLad.value;
                dmudp[offset] = muLad.derivatives[0];
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                muLad = gasPvt_.viscosity(pvtRegionIdx, TLad, pLad, RvLad);
                int offset = pu.num_phases*cellIdx + pu.phase_pos[BlackoilPhases::Vapour];
                mu[offset] = muLad.value;
                dmudp[offset] = muLad.derivatives[0];
            }
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  p      Array of n pressure values.
    /// \param[in]  T      Array of n temperature values.
    /// \param[in]  z      Array of nP surface volume values.
    /// \param[in]  cells  Array of n cell indices to be associated with the p and z values.
    /// \param[out] A      Array of nP^2 values, array must be valid before calling.
    ///                    The P^2 values for a cell give the matrix A = RB^{-1} which
    ///                    relates z to u by z = Au. The matrices are output in Fortran order.
    /// \param[out] dAdp   If non-null: array of nP^2 matrix derivative values,
    ///                    array must be valid before calling. The matrices are output
    ///                    in Fortran order.
    void BlackoilPropertiesFromDeck::matrix(const int n,
                                            const double* p,
                                            const double* T,
                                            const double* z,
                                            const int* cells,
                                            double* A,
                                            double* dAdp) const
    {
        const int np = numPhases();

        B_.resize(n*np);
        R_.resize(n*np);
        if (dAdp) {
            dB_.resize(n*np);
            dR_.resize(n*np);

            this->compute_dBdp_(n, p, T, z, cells, &B_[0], &dB_[0]);
            this->compute_dRdp_(n, p, T, z, cells, &R_[0], &dR_[0]);
        } else {
            this->compute_B_(n, p, T, z, cells, &B_[0]);
            this->compute_R_(n, p, T, z, cells, &R_[0]);
        }
        const auto& pu = phaseUsage();
        bool oil_and_gas = pu.phase_pos[BlackoilPhases::Liquid] &&
            pu.phase_pos[BlackoilPhases::Vapour];
        const int o = pu.phase_pos[BlackoilPhases::Liquid];
        const int g = pu.phase_pos[BlackoilPhases::Vapour];

        // Compute A matrix
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double* m = A + i*np*np;
            std::fill(m, m + np*np, 0.0);
            // Diagonal entries.
            for (int phase = 0; phase < np; ++phase) {
                m[phase + phase*np] = 1.0/B_[i*np + phase];
            }
            // Off-diagonal entries.
            if (oil_and_gas) {
                m[o + g*np] = R_[i*np + g]/B_[i*np + g];
                m[g + o*np] = R_[i*np + o]/B_[i*np + o];
            }
        }

        // Derivative of A matrix.
        // A     = R*inv(B) whence
        //
        // dA/dp = (dR/dp*inv(B) + R*d(inv(B))/dp)
        //       = (dR/dp*inv(B) - R*inv(B)*(dB/dp)*inv(B))
        //       = (dR/dp - A*(dB/dp)) * inv(B)
        //
        // The B matrix is diagonal and that fact is exploited in the
        // following implementation.
        if (dAdp) {
// #pragma omp parallel for
            // (1): dA/dp <- A
            std::copy(A, A + n*np*np, dAdp);

            for (int i = 0; i < n; ++i) {
                double*       m  = dAdp + i*np*np;

                // (2): dA/dp <- -dA/dp*(dB/dp) == -A*(dB/dp)
                const double* dB = & dB_[i * np];
                for (int col = 0; col < np; ++col) {
                    for (int row = 0; row < np; ++row) {
                        m[col*np + row] *= - dB[ col ]; // Note sign.
                    }
                }

                if (oil_and_gas) {
                    // (2b): dA/dp += dR/dp (== dR/dp - A*(dB/dp))
                    const double* dR = & dR_[i * np];

                    m[o*np + g] += dR[ o ];
                    m[g*np + o] += dR[ g ];
                }

                // (3): dA/dp *= inv(B) (== final result)
                const double* B = & B_[i * np];
                for (int col = 0; col < np; ++col) {
                    for (int row = 0; row < np; ++row) {
                        m[col*np + row] /= B[ col ];
                    }
                }
            }
        }
    }

    void BlackoilPropertiesFromDeck::compute_B_(const int n,
                                                const double* p,
                                                const double* T,
                                                const double* z,
                                                const int* cells,
                                                double* B) const
    {
        const auto& pu = phaseUsage();

        typedef double LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RsLad = 0.0;
        LadEval RvLad = 0.0;

        for (int i = 0; i < n; ++ i) {
            int cellIdx = cells[i];
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            pLad = p[i];
            TLad = T[i];

            int oilOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Liquid];
            int gasOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Vapour];
            int waterOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Aqua];

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                LadEval BLad = 1.0/waterPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);

                B[waterOffset] = BLad;
            }

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                double currentRs = 0.0;
                double maxRs = 0.0;
                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    currentRs = (z[oilOffset] == 0.0) ? 0.0 : z[gasOffset]/z[oilOffset];
                    maxRs = oilPvt_.saturatedGasDissolutionFactor(pvtRegionIdx, TLad, pLad);
                }
                LadEval BLad;
                if (currentRs >= maxRs) {
                    BLad = 1.0/oilPvt_.saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
                }
                else {
                    RsLad = currentRs;
                    BLad = 1.0/oilPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RsLad);
                }

                B[oilOffset] = BLad;
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                double currentRv = 0.0;
                double maxRv = 0.0;
                if (pu.phase_used[BlackoilPhases::Liquid]) {
                    currentRv = (z[gasOffset] == 0.0) ? 0.0 : z[oilOffset]/z[gasOffset];
                    maxRv = gasPvt_.saturatedOilVaporizationFactor(pvtRegionIdx, TLad, pLad);
                }
                LadEval BLad;
                if (currentRv >= maxRv) {
                    BLad = 1.0/gasPvt_.saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
                }
                else {
                    RvLad = currentRv;
                    BLad = 1.0/gasPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RvLad);
                }

                B[gasOffset] = BLad;
            }
        }
    }

    void BlackoilPropertiesFromDeck::compute_dBdp_(const int n,
                                                   const double* p,
                                                   const double* T,
                                                   const double* z,
                                                   const int* cells,
                                                   double* B,
                                                   double* dBdp) const
    {
        const auto& pu = phaseUsage();

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;
        LadEval RsLad = 0.0;
        LadEval RvLad = 0.0;

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++ i) {
            int cellIdx = cells[i];
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            pLad.value = p[i];
            TLad.value = T[i];

            int oilOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Liquid];
            int gasOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Vapour];
            int waterOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Aqua];

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                LadEval BLad = 1.0/waterPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);

                B[waterOffset] = BLad.value;
                dBdp[waterOffset] = BLad.derivatives[0];
            }

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                double currentRs = 0.0;
                double maxRs = 0.0;
                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    currentRs = (z[oilOffset] == 0.0) ? 0.0 : z[gasOffset]/z[oilOffset];
                    maxRs = oilPvt_.saturatedGasDissolutionFactor(pvtRegionIdx, TLad.value, pLad.value);
                }
                LadEval BLad;
                if (currentRs >= maxRs) {
                    BLad = 1.0/oilPvt_.saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
                }
                else {
                    RsLad.value = currentRs;
                    BLad = 1.0/oilPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RsLad);
                }

                B[oilOffset] = BLad.value;
                dBdp[oilOffset] = BLad.derivatives[0];
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                double currentRv = 0.0;
                double maxRv = 0.0;
                if (pu.phase_used[BlackoilPhases::Liquid]) {
                    currentRv = (z[gasOffset] == 0.0) ? 0.0 : z[oilOffset]/z[gasOffset];
                    maxRv = gasPvt_.saturatedOilVaporizationFactor(pvtRegionIdx, TLad.value, pLad.value);
                }
                LadEval BLad;
                if (currentRv >= maxRv) {
                    BLad = 1.0/gasPvt_.saturatedInverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad);
                }
                else {
                    RvLad.value = currentRv;
                    BLad = 1.0/gasPvt_.inverseFormationVolumeFactor(pvtRegionIdx, TLad, pLad, RvLad);
                }

                B[gasOffset] = BLad.value;
                dBdp[gasOffset] = BLad.derivatives[0];
            }
        }
    }

    void BlackoilPropertiesFromDeck::compute_R_(const int n,
                                                const double* p,
                                                const double* T,
                                                const double* z,
                                                const int* cells,
                                                double* R) const
    {
        const auto& pu = phaseUsage();

        typedef double LadEval;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;

        for (int i = 0; i < n; ++ i) {
            int cellIdx = cells[i];
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            pLad = p[i];
            TLad = T[i];

            int oilOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Liquid];
            int gasOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Vapour];
            int waterOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Aqua];

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                R[waterOffset] = 0.0; // water is always immiscible!
            }

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                LadEval RsSatLad = oilPvt_.saturatedGasDissolutionFactor(pvtRegionIdx, TLad, pLad);

                double currentRs = 0.0;
                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    currentRs = (z[oilOffset] == 0.0) ? 0.0 : z[gasOffset]/z[oilOffset];
                }

                RsSatLad = std::min(RsSatLad, currentRs);

                R[oilOffset] = RsSatLad;
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                LadEval RvSatLad = gasPvt_.saturatedOilVaporizationFactor(pvtRegionIdx, TLad, pLad);

                double currentRv = 0.0;
                if (pu.phase_used[BlackoilPhases::Liquid]) {
                    currentRv = (z[gasOffset] == 0.0) ? 0.0 : z[oilOffset]/z[gasOffset];
                }

                RvSatLad = std::min(RvSatLad, currentRv);

                R[gasOffset] = RvSatLad;
            }
        }
    }

    void BlackoilPropertiesFromDeck::compute_dRdp_(const int n,
                                                   const double* p,
                                                   const double* T,
                                                   const double* z,
                                                   const int* cells,
                                                   double* R,
                                                   double* dRdp) const
    {
        const auto& pu = phaseUsage();

        enum PressureEvalTag {};
        typedef Opm::LocalAd::Evaluation<double, PressureEvalTag, /*size=*/1> LadEval;
        typedef Opm::MathToolbox<LadEval> Toolbox;

        LadEval pLad = 0.0;
        LadEval TLad = 0.0;

        pLad.derivatives[0] = 1.0;

        for (int i = 0; i < n; ++ i) {
            int cellIdx = cells[i];
            int pvtRegionIdx = cellPvtRegionIdx_[cellIdx];
            pLad.value = p[i];
            TLad.value = T[i];

            int oilOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Liquid];
            int gasOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Vapour];
            int waterOffset = pu.num_phases*i + pu.phase_pos[BlackoilPhases::Aqua];

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                R[waterOffset] = 0.0; // water is always immiscible!
            }

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                LadEval RsSatLad = oilPvt_.saturatedGasDissolutionFactor(pvtRegionIdx, TLad, pLad);

                LadEval currentRs = 0.0;
                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    currentRs = (z[oilOffset] == 0.0) ? 0.0 : z[gasOffset]/z[oilOffset];
                }

                RsSatLad = Toolbox::min(RsSatLad, currentRs);

                R[oilOffset] = RsSatLad.value;
                dRdp[oilOffset] = RsSatLad.derivatives[0];
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                LadEval RvSatLad = gasPvt_.saturatedOilVaporizationFactor(pvtRegionIdx, TLad, pLad);

                LadEval currentRv = 0.0;
                if (pu.phase_used[BlackoilPhases::Liquid]) {
                    currentRv = (z[gasOffset] == 0.0) ? 0.0 : z[oilOffset]/z[gasOffset];
                }

                RvSatLad = Toolbox::min(RvSatLad, currentRv);

                R[gasOffset] = RvSatLad.value;
                dRdp[gasOffset] = RvSatLad.derivatives[0];
            }
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  A      Array of nP^2 values, where the P^2 values for a cell give the
    ///                    matrix A = RB^{-1} which relates z to u by z = Au. The matrices
    ///                    are assumed to be in Fortran order, and are typically the result
    ///                    of a call to the method matrix().
    /// \param[in]  cells  The index of the grid cell of each data point.
    /// \param[out] rho    Array of nP density values, array must be valid before calling.
    void BlackoilPropertiesFromDeck::density(const int n,
                                             const double* A,
                                             const int* cells,
                                             double* rho) const
    {
        const int np = numPhases();
// #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            int cellIdx = cells?cells[i]:i;
            const double *sdens = surfaceDensity(cellIdx);
            for (int phase = 0; phase < np; ++phase) {
                rho[np*i + phase] = 0.0;
                for (int comp = 0; comp < np; ++comp) {
                    rho[np*i + phase] += A[i*np*np + np*phase + comp]*sdens[comp];
                }
            }
        }
    }

    /// Densities of stock components at surface conditions.
    /// \return Array of P density values.
    const double* BlackoilPropertiesFromDeck::surfaceDensity(int cellIdx) const
    {
        const auto& pu = phaseUsage();
        int pvtRegionIdx = getTableIndex_(cellPvtRegionIndex(), cellIdx);
        return &surfaceDensities_[pvtRegionIdx*pu.num_phases];
    }

    void BlackoilPropertiesFromDeck::initSurfaceDensities_(Opm::DeckConstPtr deck)
    {
        const auto& pu = phaseUsage();
        int np = pu.num_phases;
        int numPvtRegions = 1;
        if (deck->hasKeyword("TABDIMS")) {
            const auto& tabdimsKeyword = deck->getKeyword("TABDIMS");
            numPvtRegions = tabdimsKeyword.getRecord(0).getItem("NTPVT").template get<int>(0);
        }

        const auto& densityKeyword = deck->getKeyword("DENSITY");

        surfaceDensities_.resize(np*numPvtRegions);
        for (int pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
            if (pu.phase_used[BlackoilPhases::Aqua])
                surfaceDensities_[np*pvtRegionIdx + pu.phase_pos[BlackoilPhases::Aqua]] =
                    densityKeyword.getRecord(pvtRegionIdx).getItem("WATER").getSIDouble(0);

            if (pu.phase_used[BlackoilPhases::Liquid])
                surfaceDensities_[np*pvtRegionIdx + pu.phase_pos[BlackoilPhases::Liquid]] =
                    densityKeyword.getRecord(pvtRegionIdx).getItem("OIL").getSIDouble(0);

            if (pu.phase_used[BlackoilPhases::Vapour])
                surfaceDensities_[np*pvtRegionIdx + pu.phase_pos[BlackoilPhases::Vapour]] =
                    densityKeyword.getRecord(pvtRegionIdx).getItem("GAS").getSIDouble(0);
        }
    }

    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
    /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dkr_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void BlackoilPropertiesFromDeck::relperm(const int n,
                                             const double* s,
                                             const int* cells,
                                             double* kr,
                                             double* dkrds) const
    {
        satprops_->relperm(n, s, cells, kr, dkrds);
    }


    /// \param[in]  n      Number of data points.
    /// \param[in]  s      Array of nP saturation values.
    /// \param[in]  cells  Array of n cell indices to be associated with the s values.
    /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
    /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
    ///                    array must be valid before calling.
    ///                    The P^2 derivative matrix is
    ///                           m_{ij} = \frac{dpc_i}{ds^j},
    ///                    and is output in Fortran order (m_00 m_10 m_20 m01 ...)
    void BlackoilPropertiesFromDeck::capPress(const int n,
                                              const double* s,
                                              const int* cells,
                                              double* pc,
                                              double* dpcds) const
    {
        satprops_->capPress(n, s, cells, pc, dpcds);
    }


    /// Obtain the range of allowable saturation values.
    /// In cell cells[i], saturation of phase p is allowed to be
    /// in the interval [smin[i*P + p], smax[i*P + p]].
    /// \param[in]  n      Number of data points.
    /// \param[in]  cells  Array of n cell indices.
    /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
    /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
    void BlackoilPropertiesFromDeck::satRange(const int n,
                                              const int* cells,
                                              double* smin,
                                              double* smax) const
    {
        satprops_->satRange(n, cells, smin, smax);
    }


    /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
    /// \param[in]     cell   Cell index.
    /// \param[in]     pcow   P_oil - P_water.
    /// \param[in/out] swat   Water saturation. / Possibly modified Water saturation.
    void BlackoilPropertiesFromDeck::swatInitScaling(const int cell,
                                                     const double pcow,
                                                     double & swat)
    {
        satprops_->swatInitScaling(cell, pcow, swat);
    }

} // namespace Opm

