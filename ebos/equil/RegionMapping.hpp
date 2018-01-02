// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef OPM_REGIONMAPPING_HH
#define OPM_REGIONMAPPING_HH

#include <unordered_map>
#include <vector>

namespace Ewoms {

/**
 * Forward and reverse mappings between cells and
 * regions/partitions (e.g., the ECLIPSE-style 'SATNUM',
 * 'PVTNUM', or 'EQUILNUM' arrays).
 *
 * \tparam Region Type of a forward region mapping.  Expected
 *                to provide indexed access through
 *                operator[]() as well as inner types
 *                'valueType', 'size_type', and
 *                'const_iterator'.
 */
template <class Region = std::vector<int>>
class RegionMapping
{
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

    class Range
    {
    public:
        typedef CellIter iterator;
        typedef CellIter const_iterator;

        Range()
        {};

        Range(const CellIter& beg, const CellIter& en)
            : begin_(beg)
            , end_(en)
        {};

        Range(const Range&) = default;

        CellIter& begin() { return begin_; }
        const CellIter& begin() const { return begin_; }

        const CellIter& end() const { return end_; }

        bool empty() const
        { return begin_ == end_; }

        size_t size() const
        {
            size_t ret = 0;
            for (CellIter it = begin(); it != end(); ++it)
                ++ret;

            return ret;
        }
    private:
        CellIter begin_;
        CellIter end_;
    };

    /**
     * Compute region number of given active cell.
     *
     * \param[in] c Active cell
     * \return Region to which @c c belongs.
     */
    RegionId region(const CellId c) const
    { return reg_[c]; }

    const std::vector<RegionId>&
    activeRegions() const
    {
        return rev_.active;
    }

    /**
     * Extract active cells in particular region.
     *
     * \param[in] r Region number
     *
     * \return Range of active cells in region @c r.  Empty if @c r is
     * not an active region.
     */
    Range cells(const RegionId r) const
    {
        const auto id = rev_.binid.find(r);

        if (id == rev_.binid.end()) {
            // Region 'r' not an active region.  Return empty.
            return Range(rev_.c.end(), rev_.c.end());
        }

        const auto i = id->second;

        return Range(rev_.c.begin() + rev_.p[i + 0],
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

        std::unordered_map<RegionId, Pos> binid;
        std::vector<RegionId>             active;

        std::vector<Pos>    p;   /**< Region start pointers */
        std::vector<CellId> c;   /**< Region cells */

        /**
         * Compute reverse mapping.  Standard linear insertion
         * sort algorithm.
         */
        void
        init(const Region& reg)
        {
            binid.clear();
            for (const auto& r : reg) {
                ++binid[r];
            }

            p     .clear();  p.emplace_back(0);
            active.clear();
            {
                Pos n = 0;
                for (auto& id : binid) {
                    active.push_back(id.first);
                    p     .push_back(id.second);

                    id.second = n++;
                }
            }

            for (decltype(p.size()) i = 1, sz = p.size(); i < sz; ++i) {
                p[0] += p[i];
                p[i] = p[0] - p[i];
            }

            assert (p[0] == static_cast<Pos>(reg.size()));

            c.resize(reg.size());
            {
                CellId i = 0;
                for (const auto& r : reg) {
                    auto& pos = p[binid[r] + 1];
                    c[pos++] = i++;
                }
            }

            p[0] = 0;
        }
    } rev_; /**< Reverse mapping instance */
};

} // namespace Ewoms

#endif // OPM_REGIONMAPPING_HH
