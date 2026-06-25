/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_MULTICOMM_HEADER_INCLUDED
#define OPM_MULTICOMM_HEADER_INCLUDED

#include <dune/common/hybridutilities.hh>

#if HAVE_MPI
#include <mpi.h>
#include <dune/common/parallel/mpicommunication.hh>
#endif

#include <cmath>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace Dune
{

class SeqComm
{
public:
    using size_type = std::size_t;

    static constexpr size_type size()
    {
        return 0;
    }
    template <typename T>
    void project(T& /*x*/) const
    {
    }
    template <typename T1, typename T2>
    void dot(const T1& x, const T1& y, T2& result) const
    {
        result = x.dot(y);
    }
    template <typename T>
    double norm(const T& x) const
    {
        return x.two_norm();
    }
    template <typename T>
    void copyOwnerToAll(const T& x, T& y) const
    {
        y = x;
    }

    auto communicator() const
    {
        return Dune::Communication<int>{};
    }
};

#if HAVE_MPI
class JacComm
{
public:
    using size_type = std::size_t;

    JacComm()
        : colcom_(MPI_COMM_WORLD)
    {
    }
    static constexpr size_type size()
    {
        return 0;
    }
    template <typename T>
    void project(T& /*x*/) const
    {
    }

    template <typename T1, typename T2>
    void dot(const T1& x, const T1& y, T2& result) const
    {
        result = x.dot(y);
        result = colcom_.sum(result);
    }

    template <typename T>
    double norm(const T& x) const
    {
        double result = x.dot(x);
        result = colcom_.sum(result);
        return std::sqrt(result);
    }

    template <typename T>
    void copyOwnerToAll(const T& x, T& y) const
    {
        y = x;
    }

private:
    Communication<MPI_Comm> colcom_;
};
#endif

template <typename... Args>
class MultiCommunicator : public std::tuple<Args...>
{
    using TupleType = std::tuple<Args...>;
    using type = MultiCommunicator<Args...>;
    using field_type = double;

public:
    using std::tuple<Args...>::tuple;
    using size_type = std::size_t;

    static constexpr size_type size()
    {
        return sizeof...(Args);
    }

    template <size_type index>
    typename std::tuple_element<index, TupleType>::type&
    operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable)
    {
        return std::get<index>(*this);
    }
    template <size_type index>
    const typename std::tuple_element<index, TupleType>::type&
    operator[]([[maybe_unused]] const std::integral_constant<size_type, index> indexVariable) const
    {
        return std::get<index>(*this);
    }

    template <typename T>
    void project(T& x) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) { (*this)[i].project(x[i]); });
    }
    template <typename T1, typename T2>
    void dot(const T1& x, const T1& y, T2& result) const
    {
        result = field_type(0);
        using namespace Dune::Hybrid;
        forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) {
            double result_tmp = 0;
            (*this)[i].dot(x[i], y[i], result_tmp);
            result += result_tmp;
        });
    }
    template <typename T>
    field_type norm(const T& x) const
    {
        field_type result(0);
        using namespace Dune::Hybrid;
        forEach(integralRange(Hybrid::size(*this)), [&](auto&& i) {
            double result_tmp = 0.0;
            (*this)[i].dot(x[i], x[i], result_tmp);
            result += result_tmp;
        });
        return std::sqrt(result);
    }
    template <typename T>
    void copyOwnerToAll(const T& x, T& y) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Hybrid::size(*this)),
                [&](auto&& i) { (*this)[i].copyOwnerToAll(x[i], y[i]); });
    }

    decltype(auto) communicator() const
    {
        return std::get<0>(*this).communicator();
    }
};

} // namespace Dune

#endif // OPM_MULTICOMM_HEADER_INCLUDED
