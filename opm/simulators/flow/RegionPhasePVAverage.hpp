// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor AS

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

#ifndef OPM_REGION_PHASE_POREVOL_AVERAGE_MODULE_HPP
#define OPM_REGION_PHASE_POREVOL_AVERAGE_MODULE_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace Opm {

    /// Facility for calculating volume-weighted average function values
    /// over user-defined regions in parallel.  Defaults to phase-filled
    /// pore-volume averages, but falls back to simple pore-volume averages
    /// if the saturation is identically zero in a region--i.e., calculates
    /// what the average value would be if the phase fully occupied the
    /// available pore-volume in this case.
    class RegionPhasePoreVolAverage
    {
    public:
        /// Minimal characteristics of a cell from a simulation grid.
        struct CellValue
        {
            /// Function value
            double value { 0.0 };

            /// Phase saturation
            double sat { 0.0 };

            /// Reservoir condition pore-volume
            double porv { 0.0 };
        };

        /// Compile-time disambiguation type for phases.
        struct Phase
        {
            /// Phase index
            unsigned int ix;
        };

        /// Compile-time disambiguation type for regions.
        struct Region
        {
            /// Region index
            unsigned int ix;
        };

        /// Call-back function type for accessing region arrays--typically the FIP* arrays.
        using RegionArrayAccessor = std::function<const std::vector<int>&(const std::string&)>;

        /// Constructor.
        ///
        /// \param[in] comm Grid level global communicator.
        ///
        /// \param[in] numPhases Number of phases for which to calculate
        ///   average values.
        ///
        /// \param[in] regionNames List of region sets.  Typically contains
        ///   one or more of the FIP* array names and, possibly, the PVTNUM
        ///   region set as well.
        ///
        /// \param[in] getRegionArray Call-back function for accessing
        ///   region definition from region set names.
        explicit RegionPhasePoreVolAverage(const Parallel::Communication&  comm,
                                           std::size_t                     numPhases,
                                           const std::vector<std::string>& regionNames,
                                           RegionArrayAccessor             getRegionArray);

        /// Retrieve field-level average function value for specific phase.
        ///
        /// \param[in] p Phase for which to retrieve the field-level average
        ///   function value.
        ///
        /// \return Field-level average function value for phase \p p.
        double fieldValue(const Phase& p) const;

        /// Retrieve region-level average function value for specific phase
        /// in specific region of named region set.
        ///
        /// \param[in] rset Named region set--e.g., "FIPNUM".
        ///
        /// \param[in] p Phase for which to retrieve the region-level average
        ///   function value.
        ///
        /// \param[in] r Region ID for which to retrieve the region-level
        ///   average function value.
        ///
        /// \return Region-level average function value.
        double value(std::string_view rset, const Phase& p, const Region& r) const;

        /// Clear internal arrays in preparation of accumulating
        /// region-level averages from per-cell contributions.
        void prepareAccumulation();

        /// Incorporate contributions from a single cell.
        ///
        /// \param[in] activeCell Per-rank active cell ID--typically one of
        ///   the rank's interior cells.
        ///
        /// \param[in] p Phase for which to incorporate the per-cell contribution.
        ///
        /// \param[in] cv Single cell function value contribution.
        void addCell(std::size_t activeCell, const Phase& p, const CellValue& cv);

        /// Accumulate region-level average values across MPI ranks.
        ///
        /// Typically the last step in calculating the region-level average
        /// values.  It is typically an error to call this function multiple
        /// times without an intervening call to \c prepareAccumulation().
        void accumulateParallel();

    private:
        /// Index type for value array.
        using Ix = std::vector<double>::size_type;

        /// Kinds of weighting functions.
        enum AvgType {
            SatPV,              //< Weighted by phase-filled pore-volume.
            PV,                 //< Weighted by total pore-volume (s=1).

            // ---- Must be last enumerator ----
            NumTypes,           //< Number of averge types.
        };

        /// Elements/items collected per average type.
        enum Element {
            Value,              //< Running sum of weighted function values.
            Weight,             //< Running sum of weights.

            // ---- Must be last enumerator ----
            NumElem,            //< Number of items/elements per average type.
        };

        /// MPI communication object.
        std::reference_wrapper<const Parallel::Communication> comm_;

        /// Number of phases for which to accumulate function value averages.
        std::size_t np_{};

        /// Named region sets.
        std::vector<std::string> rsetNames_{};

        /// Call-back function for accessing region set arrays.
        RegionArrayAccessor getRegionArray_;

        /// Start pointers for average values of each named region set.  In
        /// particular, \code rsStart_[i] \endcode identifies the offset
        /// into \c x_ of the first average type of the first phase of the
        /// first region in region set \code rsetNames_[i] \endcode.
        std::vector<Ix> rsStart_{};

        /// All elements--i.e., running sums--that go into calculating the
        /// field- and region-level weighted function averages per phase.
        /// We store all elements in a single linear array, and track
        /// numerators and denominators in separate elements, in order to
        /// simplify cross-rank reduction in parallel simulation runs.  See
        /// accumulateParallel() for details on this reduction process.
        ///
        /// There are \c np_ base entries for each region in each region set
        /// (+ FIELD).  Each base entry has \code AvgType::NumTypes \endcode
        /// different average types, and each average type comprises \code
        /// Element::NumElem \endcode separate elements in \c x_.
        ///
        /// You should typically access this array through the value() and
        /// weight() member functions, and the start index to both of those
        /// should be the return value from fieldStartIx() or rsetStartIx().
        std::vector<double> x_{};

        /// Compute final average value for a single region and phase.
        ///
        /// Prefers the average value weighted by phase-filled pore-volume,
        /// but falls back to average value weighted by total pore-volume if
        /// saturation is identically zero throughout the region.
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \return Weighted average function value.
        double averageValueWithFallback(Ix start) const;

        /// Compute average function value for a single region and phase.
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to compute.
        ///
        /// \return Weighted average function value.
        double averageValue(Ix start, AvgType type) const;

        /// Compute linearised value array offset for field-level average
        /// function values of a single phase.
        ///
        /// \param[in] phase Phase index for which to compute array offset.
        ///
        /// \return Value array offset for field-level averages of \p phase.
        Ix fieldStartIx(unsigned int phase) const;

        /// Compute linearised value array offset for region-level average
        /// function values of a single phase.
        ///
        /// \param[in] rset Enumerated region set.
        ///
        /// \param[in] region Region index within \p rset.
        ///
        /// \param[in] phase Phase index for which to compute array offset.
        ///
        /// \return Value array offset for region-level averages of \p phase.
        Ix rsetStartIx(std::size_t rset, int region, unsigned int phase) const;

        /// Compute linearised value array offset for average function
        /// values of a single phase.
        ///
        /// \param[in] offset Field or region base offset.
        ///
        /// \param[in] phase Phase index for which to compute array offset.
        ///
        /// \return Value array offset for weighted averages of \p phase.
        Ix startIx(std::size_t offset, unsigned int phase) const;

        /// Compute region ID of single active cell within particular region
        /// set.
        ///
        /// \param[in] rset Enumerated region set.
        ///
        /// \param[in] activeCell Per-rank active cell ID--typically one of
        ///   the rank's interior cells.
        ///
        /// \return Region ID of \p activeCell within \p rset.
        int regionIndex(std::size_t rset, std::size_t activeCell) const;

        /// Incorporate per-cell contribution into all average function types.
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] cv Single cell function value contribution.
        void add(Ix start, const CellValue& cv);

        /// Incorporate per-cell contribution into specific function type.
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to accumulate.
        ///
        /// \param[in] x Function value.
        ///
        /// \param[in] w Function weight.
        void add(Ix start, AvgType type, double x, double w);

        /// Mutable access to value item of specific average value type
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to accumulate.
        ///
        /// \return Reference to mutable element of linearised value array
        ///   corresponding to the running sum of function values.
        double& value(Ix start, AvgType type);

        /// Mutable access to weight item of specific average value type
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to accumulate.
        ///
        /// \return Reference to mutable element of linearised value array
        ///   corresponding to the running sum of function weights.
        double& weight(Ix start, AvgType type);

        /// Read-only access to value item of specific average value type
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to accumulate.
        ///
        /// \return Running sum of function value.
        double value(Ix start, AvgType type) const;

        /// Read-only access to weight item of specific average value type
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to accumulate.
        ///
        /// \return Running sum of function weights.
        double weight(Ix start, AvgType type) const;

        /// Compute value array index of particular item of particular
        /// average value
        ///
        /// \param[in] start Offset into linearised value array (\c x_)
        ///   corresponding to the first average value type of a particular
        ///   phase in a particular region.  Usually calculated by
        ///   fieldStartIx() or rsetStartIx().
        ///
        /// \param[in] type Which kind of average value to reference.
        ///
        /// \param[in] element Which running sum element to reference for
        ///   this particular average value.
        ///
        /// \return Index into linearised value array \c x_ corresponding to
        ///   this particular element of this particular average function
        ///   type.
        Ix valueArrayIndex(Ix start, AvgType type, Element element) const;
    };
} // namespace Opm

#endif // OPM_REGION_PHASE_POREVOL_AVERAGE_MODULE_HPP
