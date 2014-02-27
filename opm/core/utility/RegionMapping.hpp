/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_REGIONMAPPING_HEADER_INCLUDED
#define OPM_REGIONMAPPING_HEADER_INCLUDED

#include <vector>

namespace Opm
{

    /**
     * Forward and reverse mappings between cells and
     * regions/partitions (e.g., the ECLIPSE-style 'SATNUM',
     * 'PVTNUM', or 'EQUILNUM' arrays).
     *
     * \tparam Region Type of a forward region mapping.  Expected
     *                to provide indexed access through
     *                operator[]() as well as inner types
     *                'value_type', 'size_type', and
     *                'const_iterator'.
     */
    template < class Region = std::vector<int> >
    class RegionMapping {
    public:
        /**
         * Constructor.
         *
         * \param[in] reg Forward region mapping, restricted to
         *                active cells only.
         */
        explicit
        RegionMapping(const Region& reg)
        : reg_(reg)
        {
            rev_.init(reg_);
        }

        /**
         * Type of forward (cell-to-region) mapping result.
         * Expected to be an integer.
         */
        typedef typename Region::value_type RegionId;

        /**
         * Type of reverse (region-to-cell) mapping (element)
         * result.
         */
        typedef typename Region::size_type CellId;

        /**
         * Type of reverse region-to-cell range bounds and
         * iterators.
         */
        typedef typename std::vector<CellId>::const_iterator CellIter;

        /**
         * Range of cells.  Result from reverse (region-to-cell)
         * mapping.
         */
        class CellRange {
        public:
            /**
             * Constructor.
             *
             * \param[in] b Beginning of range.
             * \param[in] e One past end of range.
             */
            CellRange(const CellIter b,
                      const CellIter e)
                : b_(b), e_(e)
            {}

            /**
             * Read-only iterator on cell ranges.
             */
            typedef CellIter const_iterator;

            /**
             * Size type for this range.
             */
            typedef typename std::vector<CellId>::size_type size_type;

            /**
             * Beginning of cell range.
             */
            const_iterator begin() const { return b_; }

            /**
             * One past end of cell range.
             */
            const_iterator end()   const { return e_; }

            /**
             * Number of elements in the range.
             */
            size_type size() const { return e_ - b_; }

        private:
            const_iterator b_;
            const_iterator e_;
        };

        /**
         * Number of declared regions in cell-to-region mapping.
         */
        RegionId
        numRegions() const { return RegionId(rev_.p.size()) - 1; }

        /**
         * Compute region number of given active cell.
         *
         * \param[in] c Active cell
         * \return Region to which @c c belongs.
         */
        RegionId
        region(const CellId c) const { return reg_[c]; }

        /**
         * Extract active cells in particular region.
         *
         * \param[in] r Region number
         * \returns Range of active cells in region @c r.
         */
        CellRange
        cells(const RegionId r) const {
            const RegionId i = r - rev_.low;
            return CellRange(rev_.c.begin() + rev_.p[i + 0],
                             rev_.c.begin() + rev_.p[i + 1]);
        }

    private:
        /**
         * Copy of forward region mapping (cell-to-region).
         */
        Region reg_;

        /**
         * Reverse mapping (region-to-cell).
         */
        struct {
            typedef typename std::vector<CellId>::size_type Pos;
            std::vector<Pos>    p;   /**< Region start pointers */
            std::vector<CellId> c;   /**< Region cells */
            RegionId            low; /**< Smallest region number */

            /**
             * Compute reverse mapping.  Standard linear insertion
             * sort algorithm.
             */
            void
            init(const Region& reg)
            {
                typedef typename Region::const_iterator CI;
                const std::pair<CI,CI>
                    m = std::minmax_element(reg.begin(), reg.end());

                low  = *m.first;

                const typename Region::size_type
                    n = *m.second - low + 1;

                p.resize(n + 1);  std::fill(p.begin(), p.end(), Pos(0));
                for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                    p[ reg[i] - low + 1 ] += 1;
                }

                for (typename std::vector<Pos>::size_type
                         i = 1, sz = p.size(); i < sz; ++i) {
                    p[0] += p[i];
                    p[i]  = p[0] - p[i];
                }

                assert (p[0] ==
                        static_cast<typename Region::size_type>(reg.size()));

                c.resize(reg.size());
                for (CellId i = 0, nc = reg.size(); i < nc; ++i) {
                    c[ p[ reg[i] - low + 1 ] ++ ] = i;
                }

                p[0] = 0;
            }
        } rev_; /**< Reverse mapping instance */
    };

} // namespace Opm

#endif // OPM_REGIONMAPPING_HEADER_INCLUDED
