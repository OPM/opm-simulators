/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUISTL_DETAIL_HAS_FUNCTION_HPP
#define OPM_CUISTL_DETAIL_HAS_FUNCTION_HPP
#include <type_traits>
/**
 * Simple utility structs to test for the existence of functions in types.
 *
 * Note that there are alternatives to this, see for instance
 * https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function , however,
 * this is by far the cleanest approach for where this is going to be used for now.
 *
 * TODO: Use the requires-keyword once C++20 becomes availble (https://en.cppreference.com/w/cpp/language/constraints ).
 * With C++20 this file can be removed.
 */
namespace Opm::cuistl::detail
{

/**
 * @brief The has_should_call_pre class detects the presence of the method shouldCallPre
 *
 * Usage:
 *
 * @code{.cpp}
 * if constexpr (has_should_call_pre<decltype(preconditioner)>::value) {
 *     // We know that the function shouldCallPre is present:
 *     auto shouldCallPre = preconditioner.shouldCallPre();
 * }
 * @endcode
 *
 * @note This is mainly done in the GPU preconditioner to avoid having to copy data in the pre step.
 */
template <typename T>
class has_should_call_pre
{
    template <typename U>
    static std::true_type test(decltype(&U::shouldCallPre));
    template <typename U>
    static std::false_type test(...);

public:
    static constexpr bool value = std::is_same_v<decltype(test<T>(0)), std::true_type>;
};

/**
 * @brief The has_should_call_post class detects the presence of the method shouldCallPost
 *
 * Usage:
 *
 * @code{.cpp}
 * if constexpr (has_should_call_post<decltype(preconditioner)>::value) {
 *     // We know that the function shouldCallPost is present:
 *     auto shouldCallPost = preconditioner.shouldCallPost();
 * }
 * @endcode
 *
 * @note This is mainly done in the GPU preconditioner to avoid having to copy data in the post step.
 */
template <typename T>
class has_should_call_post
{
    template <typename U>
    static std::true_type test(decltype(&U::shouldCallPost));
    template <typename U>
    static std::false_type test(...);

public:
    static constexpr bool value = std::is_same_v<decltype(test<T>(0)), std::true_type>;
};


/**
 * @brief The has_communication class checks if the type has the member function getCommunication.
 *
 * This is used in the SolverAdapter to check if the operator has a communication, and it will then select a different
 * path accordingly.
 */
template <typename T>
class has_communication
{
    template <typename U>
    static std::true_type test(decltype(&U::getCommunication));
    template <typename U>
    static std::false_type test(...);

public:
    static constexpr bool value = std::is_same_v<decltype(test<T>(0)), std::true_type>;
};

/**
 * @brief The is_a_well_operator class tries to guess if the operator is a well type operator
 *
 * @note This is mainly done in the solver adapter to detect incompatible runtime arguments. When the GPU linear solve
 * paths supports wells, this class can be removed.
 */
template <typename T>
class is_a_well_operator
{
    template <typename U>
    static std::true_type test(decltype(&U::addWellPressureEquations));
    template <typename U>
    static std::false_type test(...);

public:
    static constexpr bool value = std::is_same_v<decltype(test<T>(0)), std::true_type>;
};

} // namespace Opm::cuistl::detail
#endif
