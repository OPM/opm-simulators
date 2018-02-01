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
 *
 * \brief This is a simple test that illustrates how to use the Opm::ConditionalStorage
 *        class.
*/
#include "config.h"

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <string>
#include <cstdlib>
#include <type_traits>

template <bool foo>
class EnsureCompileTimeConstant
{};

class IAmAnIslandLeaveMeAlone
{
public:
    IAmAnIslandLeaveMeAlone(int, int)
    {}

    IAmAnIslandLeaveMeAlone& operator=(const IAmAnIslandLeaveMeAlone&) = delete;
};

int main()
{
    {
        typedef Opm::ConditionalStorage<true, std::string> ConditionalTrueString;
        ConditionalTrueString foo; // default constructor
        ConditionalTrueString bar("hello"); // construct using arguments
        ConditionalTrueString baz(bar); // copy constructor

        EnsureCompileTimeConstant<ConditionalTrueString::condition> hello OPM_UNUSED;
        if (!std::is_same<typename ConditionalTrueString::type, std::string>::value)
            // something went wrong with the exported type
            std::abort();

        if (ConditionalTrueString::condition != true)
            // the condition is not exported correctly
            std::abort();

        if (*bar != "hello")
            // value constructor did not work
            std::abort();

        if (*bar != *baz)
            // copy constructor did not work
            std::abort();

        if (foo->size() != 0)
            // default constructor did not work
            std::abort();

        // the assignment operator for the "wrapper" object should work
        foo = baz;
        if (*foo != *baz)
            // assignment operator did not work
            std::abort();
    }

    {
        typedef Opm::ConditionalStorage<false, std::string> ConditionalFalseString;
        ConditionalFalseString foo; // default constructor
        ConditionalFalseString bar("hello"); // construct by value
        ConditionalFalseString OPM_UNUSED baz(bar); // copy constructor

        EnsureCompileTimeConstant<ConditionalFalseString::condition> hello OPM_UNUSED;
        if (!std::is_same<typename ConditionalFalseString::type, std::string>::value)
            // something went wrong with the exported type
            std::abort();

        if (ConditionalFalseString::condition != false)
            // the condition is not exported correctly
            std::abort();

        // the assignment operator for the "wrapper" object should always work
        baz = foo;

        try {
            *bar;

            // this is supposed to throw an std::logic_error
            std::abort();
        }
        catch (std::logic_error &) {}

        try {
            foo->size();

            // this is supposed to throw an std::logic_error
            std::abort();
        }
        catch (std::logic_error &) {}
    }

    {
        typedef Opm::ConditionalStorage<true, IAmAnIslandLeaveMeAlone> ConditionalTrueIsland;
        ConditionalTrueIsland OPM_UNUSED foo(1, 2);
        // ConditionalTrueIsland OPM_UNUSED bar; // compiler fails because of missing default ctor
    }

    {
        typedef Opm::ConditionalStorage<false, IAmAnIslandLeaveMeAlone> ConditionalFalseIsland;
        ConditionalFalseIsland OPM_UNUSED foo(1, 2);
        // ConditionalFalseIsland OPM_UNUSED bar; // compiler fails because of missing default ctor
    }

    return 0;
}
