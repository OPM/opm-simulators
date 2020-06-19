// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
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
/*!
 * \file
 * \copydoc Opm::ThreadedEntityIterator
 */
#ifndef EWOMS_THREADED_ENTITY_ITERATOR_HH
#define EWOMS_THREADED_ENTITY_ITERATOR_HH

#include <thread>
#include <mutex>

namespace Opm {

/*!
 * \brief Provides an STL-iterator like interface to iterate over the enties of a
 *        GridView in OpenMP threaded applications
 *
 * ATTENTION: This class must be instantiated in a sequential context!
 */
template <class GridView, int codim>
class ThreadedEntityIterator
{
    using Entity = typename GridView::template Codim<codim>::Entity;
    using EntityIterator = typename GridView::template Codim<codim>::Iterator;
public:
    ThreadedEntityIterator(const GridView& gridView)
        : gridView_(gridView)
        , sequentialIt_(gridView_.template begin<codim>())
        , sequentialEnd_(gridView.template end<codim>())
    { }

    ThreadedEntityIterator(const ThreadedEntityIterator& other) = default;

    // begin iterating over the grid in parallel
    EntityIterator beginParallel()
    {
        mutex_.lock();
        auto tmp = sequentialIt_;
        if (sequentialIt_ != sequentialEnd_)
            ++sequentialIt_; // make the next thread look at the next element
        mutex_.unlock();

        return tmp;
    }

    // returns true if the last element was reached
    bool isFinished(const EntityIterator& it) const
    { return it == sequentialEnd_; }

    // make sure that the loop over the grid is finished
    void setFinished()
    {
        mutex_.lock();
        sequentialIt_ = sequentialEnd_;
        mutex_.unlock();
    }

    // prefix increment: goes to the next element which is not yet worked on by any
    // thread
    EntityIterator increment()
    {
        mutex_.lock();
        auto tmp = sequentialIt_;
        if (sequentialIt_ != sequentialEnd_)
            ++sequentialIt_;
        mutex_.unlock();

        return tmp;
    }

private:
    GridView gridView_;
    EntityIterator sequentialIt_;
    EntityIterator sequentialEnd_;

    std::mutex mutex_;
};
} // namespace Opm

#endif
