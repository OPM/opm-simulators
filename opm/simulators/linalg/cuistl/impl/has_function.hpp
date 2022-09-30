/*
  Copyright SINTEF AS 2022

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
#ifndef OPM_CUISTL_IMPL_HAS_FUNCTION_HEADER
#define OPM_CUISTL_IMPL_HAS_FUNCTION_HEADER
namespace Opm::cuistl::impl
{
template <typename T>
class has_communication
{
    using yes_type = char;
    using no_type = long;
    template <typename U>
    static yes_type test(decltype(&U::getCommunication));
    template <typename U>
    static no_type test(...);

public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);
};

template <typename T>
class is_a_well_operator
{
    using yes_type = char;
    using no_type = long;
    template <typename U>
    static yes_type test(decltype(&U::addWellPressureEquations));
    template <typename U>
    static no_type test(...);

public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);
};

template <typename T>
class has_should_call_pre
{
    using yes_type = char;
    using no_type = long;
    template <typename U>
    static yes_type test(decltype(&U::shouldCallPre));
    template <typename U>
    static no_type test(...);

public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);
};

template <typename T>
class has_should_call_post
{
    using yes_type = char;
    using no_type = long;
    template <typename U>
    static yes_type test(decltype(&U::shouldCallPost));
    template <typename U>
    static no_type test(...);

public:
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);
};

} // namespace Opm::cuistl::impl
#endif
