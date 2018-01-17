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

#ifndef OPM_BLACKOILPROPERTIESFROMDECK_HEADER_INCLUDED
#define OPM_BLACKOILPROPERTIESFROMDECK_HEADER_INCLUDED


#include <opm/core/props/BlackoilPropertiesInterface.hpp>
#include <opm/core/props/rock/RockFromDeck.hpp>
#include <opm/core/props/satfunc/SaturationPropsFromDeck.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <memory>

struct UnstructuredGrid;

namespace Opm
{

    /// Concrete class implementing the blackoil property interface,
    /// reading all data and properties from eclipse deck input.
    class BlackoilPropertiesFromDeck : public BlackoilPropertiesInterface
    {
    public:
        typedef typename SaturationPropsFromDeck::MaterialLawManager MaterialLawManager;

        /// Initialize from deck and grid.
        /// \param[in]  deck     Deck input parser
        /// \param[in]  grid     Grid to which property object applies, needed for the
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        BlackoilPropertiesFromDeck(const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   const UnstructuredGrid& grid, bool init_rock=true );

        /// Initialize from deck, grid and parameters.
        /// \param[in]  deck     Deck input parser
        /// \param[in]  grid     Grid to which property object applies, needed for the
        ///                      mapping from cell indices (typically from a processed grid)
        ///                      to logical cartesian indices consistent with the deck.
        /// \param[in]  param    Parameters. Accepted parameters include:
        ///                        pvt_tab_size (200)          number of uniform sample points for dead-oil pvt tables.
        ///                        sat_tab_size (200)          number of uniform sample points for saturation tables.
        ///                        threephase_model("simple")  three-phase relperm model (accepts "simple" and "stone2").
        ///                      For both size parameters, a 0 or negative value indicates that no spline fitting is to
        ///                      be done, and the input fluid data used directly for linear interpolation.
        BlackoilPropertiesFromDeck(const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   const UnstructuredGrid& grid,
                                   const ParameterGroup& param,
                                   bool init_rock=true);

        BlackoilPropertiesFromDeck(const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   int number_of_cells,
                                   const int* global_cell,
                                   const int* cart_dims,
                                   bool init_rock=true);

        BlackoilPropertiesFromDeck(const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   int number_of_cells,
                                   const int* global_cell,
                                   const int* cart_dims,
                                   const ParameterGroup& param,
                                   bool init_rock=true);

        BlackoilPropertiesFromDeck(const Opm::Deck& deck,
                                   const Opm::EclipseState& eclState,
                                   std::shared_ptr<MaterialLawManager> materialLawManager,
                                   int number_of_cells,
                                   const int* global_cell,
                                   const int* cart_dims,
                                   const ParameterGroup& param,
                                   bool init_rock=true);

        /// Destructor.
        virtual ~BlackoilPropertiesFromDeck();


        // ---- Rock interface ----

        /// \return   D, the number of spatial dimensions.
        virtual int numDimensions() const;

        /// \return   N, the number of cells.
        virtual int numCells() const;

        /// Return an array containing the PVT table index for each
        /// grid cell
        virtual const int* cellPvtRegionIndex() const
        { return &cellPvtRegionIdx_[0]; }

        /// \return   Array of N porosity values.
        virtual const double* porosity() const;

        /// \return   Array of ND^2 permeability values.
        ///           The D^2 permeability values for a cell are organized as a matrix,
        ///           which is symmetric (so ordering does not matter).
        virtual const double* permeability() const;


        // ---- Fluid interface ----

        /// \return   P, the number of phases (also the number of components).
        virtual int numPhases() const;

        /// \return   Object describing the active phases.
        virtual PhaseUsage phaseUsage() const;

        /// \param[in]  n      Number of data points.
        /// \param[in]  p      Array of n pressure values.
        /// \param[in]  T      Array of n temperature values.
        /// \param[in]  z      Array of nP surface volume values.
        /// \param[in]  cells  Array of n cell indices to be associated with the p and z values.
        /// \param[out] mu     Array of nP viscosity values, array must be valid before calling.
        /// \param[out] dmudp  If non-null: array of nP viscosity derivative values,
        ///                    array must be valid before calling.
        virtual void viscosity(const int n,
                               const double* p,
                               const double* T,
                               const double* z,
                               const int* cells,
                               double* mu,
                               double* dmudp) const;

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
        virtual void matrix(const int n,
                            const double* p,
                            const double* T,
                            const double* z,
                            const int* cells,
                            double* A,
                            double* dAdp) const;


        /// Densities of stock components at reservoir conditions.
        /// \param[in]  n      Number of data points.
        /// \param[in]  A      Array of nP^2 values, where the P^2 values for a cell give the
        ///                    matrix A = RB^{-1} which relates z to u by z = Au. The matrices
        ///                    are assumed to be in Fortran order, and are typically the result
        ///                    of a call to the method matrix().
        /// \param[in]  cells  The index of the grid cell of each data point.
        /// \param[out] rho    Array of nP density values, array must be valid before calling.
        virtual void density(const int n,
                             const double* A,
                             const int* cells,
                             double* rho) const;

        /// Densities of stock components at surface conditions.
        /// \return Array of P density values.
        virtual const double* surfaceDensity(int cellIdx = 0) const;

        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] kr     Array of nP relperm values, array must be valid before calling.
        /// \param[out] dkrds  If non-null: array of nP^2 relperm derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dkr_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)
        virtual void relperm(const int n,
                             const double* s,
                             const int* cells,
                             double* kr,
                             double* dkrds) const;


        /// \param[in]  n      Number of data points.
        /// \param[in]  s      Array of nP saturation values.
        /// \param[in]  cells  Array of n cell indices to be associated with the s values.
        /// \param[out] pc     Array of nP capillary pressure values, array must be valid before calling.
        /// \param[out] dpcds  If non-null: array of nP^2 derivative values,
        ///                    array must be valid before calling.
        ///                    The P^2 derivative matrix is
        ///                           m_{ij} = \frac{dpc_i}{ds^j},
        ///                    and is output in Fortran order (m_00 m_10 m_20 m_01 ...)
        virtual void capPress(const int n,
                              const double* s,
                              const int* cells,
                              double* pc,
                              double* dpcds) const;


        /// Obtain the range of allowable saturation values.
        /// In cell cells[i], saturation of phase p is allowed to be
        /// in the interval [smin[i*P + p], smax[i*P + p]].
        /// \param[in]  n      Number of data points.
        /// \param[in]  cells  Array of n cell indices.
        /// \param[out] smin   Array of nP minimum s values, array must be valid before calling.
        /// \param[out] smax   Array of nP maximum s values, array must be valid before calling.
        virtual void satRange(const int n,
                              const int* cells,
                              double* smin,
                              double* smax) const;
                      
                     
        /// Update capillary pressure scaling according to pressure diff. and initial water saturation.
        /// \param[in]     cell   Cell index. 
        /// \param[in]     pcow   P_oil - P_water.
        /// \param[in/out] swat   Water saturation. / Possibly modified Water saturation.
        virtual  void swatInitScaling(const int cell,
                                      const double pcow, 
                                      double & swat);

        const OilPvtMultiplexer<double>& oilPvt() const
        {
            return oilPvt_;
        }

        const GasPvtMultiplexer<double>& gasPvt() const
        {
            return gasPvt_;
        }

        const WaterPvtMultiplexer<double>& waterPvt() const
        {
            return waterPvt_;
        }

    private:
        int getTableIndex_(const int* pvtTableIdx, int cellIdx) const
        {
            if (!pvtTableIdx)
                return 0;
            return pvtTableIdx[cellIdx];
        }

        void initSurfaceDensities_(const Opm::Deck& deck);

        void compute_B_(const int n,
                        const double* p,
                        const double* T,
                        const double* z,
                        const int* cells,
                        double* B) const;

        void compute_dBdp_(const int n,
                           const double* p,
                           const double* T,
                           const double* z,
                           const int* cells,
                           double* B,
                           double* dBdp) const;

        void compute_R_(const int n,
                        const double* p,
                        const double* T,
                        const double* z,
                        const int* cells,
                        double* R) const;

        void compute_dRdp_(const int n,
                           const double* p,
                           const double* T,
                           const double* z,
                           const int* cells,
                           double* R,
                           double* dRdp) const;

        void init(const Opm::Deck& deck,
                  const Opm::EclipseState& eclState,
                  std::shared_ptr<MaterialLawManager> materialLawManager,
                  int number_of_cells,
                  const int* global_cell,
                  const int* cart_dims,
                  bool init_rock);

        void init(const Opm::Deck& deck,
                  const Opm::EclipseState& eclState,
                  std::shared_ptr<MaterialLawManager> materialLawManager,
                  int number_of_cells,
                  const int* global_cell,
                  const int* cart_dims,
                  const ParameterGroup& param,
                  bool init_rock);

        RockFromDeck rock_;
        PhaseUsage phaseUsage_;
        std::vector<int> cellPvtRegionIdx_;
        OilPvtMultiplexer<double> oilPvt_;
        GasPvtMultiplexer<double> gasPvt_;
        WaterPvtMultiplexer<double> waterPvt_;
        std::shared_ptr<MaterialLawManager> materialLawManager_;
        std::shared_ptr<SaturationPropsInterface> satprops_;
        std::vector<double> surfaceDensities_;
        mutable std::vector<double> B_;
        mutable std::vector<double> dB_;
        mutable std::vector<double> R_;
        mutable std::vector<double> dR_;
    };



} // namespace Opm

#endif // OPM_BLACKOILPROPERTIESFROMDECK_HEADER_INCLUDED
