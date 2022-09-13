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
#ifndef ECL_MPI_SERIALIZER_HH
#define ECL_MPI_SERIALIZER_HH

#include <opm/simulators/utils/MPIPacker.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <map>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

namespace Opm {

namespace detail
{

template<typename ...Ts>
struct MakeVariantImpl
{

template<std::size_t Index, typename, typename ...Rest>
static decltype(auto) make_variant(std::size_t index)
{
  if(Index == index)
    return std::variant<Ts...>{std::in_place_index_t<Index>{}};

  if constexpr(sizeof...(Rest) != 0)
    return make_variant<Index + 1, Rest...>(index);
  else
    throw std::runtime_error("Invalid variant index");
}

};

template<typename ...Ts>
decltype(auto) make_variant(std::size_t index)
{
  return detail::MakeVariantImpl<Ts...>::template make_variant<0, Ts...>(index);
}

template<class T>
using remove_cvr_t = std::remove_cv_t<std::remove_reference_t<T>>;

} // namespace detail

/*! \brief Class for (de-)serializing and broadcasting data in parallel.
 *!  \details If the class has a serializeOp member this is used,
 *            if not it is passed on to the underlying primitive serializer.
*/

class EclMpiSerializer {
public:
    //! \brief Constructor.
    //! \param comm The global communicator to broadcast using
    explicit EclMpiSerializer(Opm::Parallel::Communication comm) :
        m_comm(comm)
    {}

    //! \brief (De-)serialization for simple types.
    //! \details The data handled by this depends on the underlying serialization used.
    //!          Currently you can call this for scalars, and stl containers with scalars.
    template<class T>
    void operator()(const T& data)
    {
        if constexpr (is_ptr<T>::value) {
            ptr(data);
        } else if constexpr (is_pair<T>::value) {
            pair(data);
        } else if constexpr (is_variant<T>::value) {
            variant(data);
        } else if constexpr (is_optional<T>::value) {
            optional(data);
        } else if constexpr (is_vector<T>::value) {
            vector(const_cast<T&>(data));
        } else if constexpr (is_map<T>::value) {
            map(const_cast<T&>(data));
        } else if constexpr (is_array<T>::value) {
            array(const_cast<T&>(data));
        } else if constexpr (is_set<T>::value) {
            set(const_cast<T&>(data));
        } else if constexpr (has_serializeOp<detail::remove_cvr_t<T>>::value) {
          const_cast<T&>(data).serializeOp(*this);
        } else {
          if (m_op == Operation::PACKSIZE)
              m_packSize += Mpi::packSize(data, m_comm);
          else if (m_op == Operation::PACK)
              Mpi::pack(data, m_buffer, m_position, m_comm);
          else if (m_op == Operation::UNPACK)
              Mpi::unpack(const_cast<T&>(data), m_buffer, m_position, m_comm);
        }
    }

    //! \brief Handler for vectors.
    //! \tparam T Type for vector elements
    //! \param data The vector to (de-)serialize
    template <typename T>
    void vector(std::vector<T>& data)
    {
        auto handle = [&](auto& d)
        {
            for (auto& it : d) {
              if constexpr (is_pair<T>::value)
                  pair(it);
              else if constexpr (is_ptr<T>::value)
                  ptr(it);
              else if constexpr (is_vector<T>::value)
                vector(it);
              else if constexpr (has_serializeOp<T>::value)
                  it.serializeOp(*this);
              else
                  (*this)(it);
            }
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            handle(data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            handle(data);
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            data.resize(size);
            handle(data);
        }
    }

    //! \brief Handler for bool vectors.
    //! \param data The vector to (de-)serialize
    void vector(std::vector<bool>& data)
    {
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            m_packSize += data.size()*Mpi::packSize(bool(), m_comm);
        } else if (m_op == Operation::PACK) {
            (*this)(data.size());
            for (const auto entry : data) { // Not a reference: vector<bool> range
                bool b = entry;
                (*this)(b);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            (*this)(size);
            data.clear();
            data.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                bool entry;
                (*this)(entry);
                data.push_back(entry);
            }
        }
    }

    template <class Array>
    void array(Array& data)
    {
        using T = typename Array::value_type;

        auto handle = [&](auto& d) {
            for (auto& it : d) {
                if constexpr (is_pair<T>::value)
                    pair(it);
                else if constexpr (is_ptr<T>::value)
                    ptr(it);
                else if constexpr (has_serializeOp<T>::value)
                    it.serializeOp(*this);
                else
                    (*this)(it);
            }
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            handle(data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            handle(data);
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            handle(data);
        }
    }

    //! \brief Handler for std::variant.
    //! \param data The variant to (de-)serialize
    template<class... Args>
    void variant(const std::variant<Args...>& data)
    {
        auto visitor = [&](auto& d)
        {
            if constexpr (has_serializeOp<detail::remove_cvr_t<decltype(d)>>::value)
                const_cast<detail::remove_cvr_t<decltype(d)>&>(d).serializeOp(*this);
            else
                (*this)(d);
        };
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.index(), m_comm);
            std::visit(visitor, data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.index(), m_buffer, m_position, m_comm);
            std::visit(visitor, data);
        } else if (m_op == Operation::UNPACK) {
            size_t index;
            Mpi::unpack(index, m_buffer, m_position, m_comm);
            auto& data_mut = const_cast<std::variant<Args...>&>(data);
            data_mut = detail::make_variant<Args...>(index);
            std::visit(visitor, data_mut);
        }
    }

    //! \brief Handler for std::optional.
    //! \tparam T Type for data
    //! \param data The optional to (de-)serialize
    template<class T>
    void optional(const std::optional<T>& data)
    {
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.has_value(), m_comm);
            if (data.has_value()) {
                (*this)(*data);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.has_value(), m_buffer, m_position, m_comm);
            if (data.has_value()) {
                (*this)(*data);
            }
        } else if (m_op == Operation::UNPACK) {
            bool has;
            Mpi::unpack(has, m_buffer, m_position, m_comm);
            if (has) {
                T res;
                (*this)(res);
                const_cast<std::optional<T>&>(data) = res;
            }
        }

    }

    //! \brief Handler for maps.
    //! \tparam Map map type
    //! \param map The map to (de-)serialize
    template<class Map>
    void map(Map& data)
    {
        using Key = typename Map::key_type;
        using Data = typename Map::mapped_type;

        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                vector(d);
            else if constexpr (is_ptr<Data>::value)
                ptr(d);
            else if constexpr (is_map<Data>::value)
                map(d);
            else if constexpr (has_serializeOp<Data>::value)
                d.serializeOp(*this);
            else
                (*this)(d);
        };

        auto keyHandle = [&](auto& d)
        {
              if constexpr (is_pair<Key>::value)
                  pair(d);
              else if constexpr (has_serializeOp<Key>::value)
                  d.serializeOp(*this);
              else
                  (*this)(d);
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            for (auto& it : data) {
                keyHandle(it.first);
                handle(it.second);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            for (auto& it : data) {
                keyHandle(it.first);
                handle(it.second);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            for (size_t i = 0; i < size; ++i) {
                Key key;
                keyHandle(key);
                Data entry;
                handle(entry);
                data.insert(std::make_pair(key, entry));
            }
        }
    }

    template<class Set>
    void set(Set& data)
    {
        using Data = typename Set::value_type;

        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                vector(d);
            else if constexpr (is_ptr<Data>::value)
                ptr(d);
            else if constexpr (has_serializeOp<Data>::value)
                d.serializeOp(*this);
            else
                (*this)(d);
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), m_comm);
            for (auto& it : data) {
                handle(it);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, m_comm);
            for (auto& it : data) {
                handle(it);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, m_comm);
            for (size_t i = 0; i < size; ++i) {
                Data entry;
                handle(entry);
                data.insert(entry);
            }
        }
    }

    //! \brief Call this to serialize data.
    //! \tparam T Type of class to serialize
    //! \param data Class to serialize
    template<class T>
    void pack(T& data)
    {
        m_op = Operation::PACKSIZE;
        m_packSize = 0;
        data.serializeOp(*this);
        m_position = 0;
        m_buffer.resize(m_packSize);
        m_op = Operation::PACK;
        data.serializeOp(*this);
    }

    //! \brief Call this to de-serialize data.
    //! \tparam T Type of class to de-serialize
    //! \param data Class to de-serialize
    template<class T>
    void unpack(T& data)
    {
        m_position = 0;
        m_op = Operation::UNPACK;
        data.serializeOp(*this);
    }

    //! \brief Serialize and broadcast on root process, de-serialize on
    //! others.
    //!
    //! \tparam T Type of class to broadcast
    //! \param data Class to broadcast
    //! \param root Process to broadcast from
    template<class T>
    void broadcast(T& data, int root = 0)
    {
        if (m_comm.size() == 1)
            return;

        if (m_comm.rank() == root) {
            try {
                pack(data);
                m_packSize = m_position;
                m_comm.broadcast(&m_packSize, 1, root);
                m_comm.broadcast(m_buffer.data(), m_position, root);
            } catch (...) {
                m_packSize = std::numeric_limits<size_t>::max();
                m_comm.broadcast(&m_packSize, 1, root);
                throw;
            }
        } else {
            m_comm.broadcast(&m_packSize, 1, root);
            if (m_packSize == std::numeric_limits<size_t>::max()) {
                throw std::runtime_error("Error detected in parallel serialization");
            }

            m_buffer.resize(m_packSize);
            m_comm.broadcast(m_buffer.data(), m_packSize, root);
            unpack(data);
        }
    }

    template<typename... Args>
    void broadcast(int root, Args&&... args)
    {
        if (m_comm.size() == 1)
            return;

        if (m_comm.rank() == root) {
            try {
                m_op = Operation::PACKSIZE;
                m_packSize = 0;
                variadic_call(args...);
                m_position = 0;
                m_buffer.resize(m_packSize);
                m_op = Operation::PACK;
                variadic_call(args...);
                m_packSize = m_position;
                m_comm.broadcast(&m_packSize, 1, 0);
                m_comm.broadcast(m_buffer.data(), m_position, 0);
            } catch (...) {
                m_packSize = std::numeric_limits<size_t>::max();
                m_comm.broadcast(&m_packSize, 1, 0);
                throw;
            }
        } else {
            m_comm.broadcast(&m_packSize, 1, 0);
            if (m_packSize == std::numeric_limits<size_t>::max()) {
                throw std::runtime_error("Error detected in parallel serialization");
            }
            m_buffer.resize(m_packSize);
            m_comm.broadcast(m_buffer.data(), m_packSize, 0);
            m_position = 0;
            m_op = Operation::UNPACK;
            variadic_call(std::forward<Args>(args)...);
        }
    }

    //! \brief Serialize and broadcast on root process, de-serialize and append on
    //! others.
    //!
    //! \tparam T Type of class to broadcast
    //! \param data Class to broadcast
    //! \param root Process to broadcast from
    template<class T>
    void append(T& data, int root = 0)
    {
        if (m_comm.size() == 1)
            return;

        T tmp;
        T& bcast = m_comm.rank() == root ? data : tmp;
        broadcast(bcast);

        if (m_comm.rank() != root)
            data.append(tmp);
    }

    //! \brief Returns current position in buffer.
    size_t position() const
    {
        return m_position;
    }

    //! \brief Returns true if we are currently doing a serialization operation.
    bool isSerializing() const
    {
        return m_op != Operation::UNPACK;
    }

protected:
    template<typename T, typename... Args>
    void variadic_call(T& first,
                       Args&&... args)
    {
      (*this)(first);
      if constexpr (sizeof...(args) > 0)
          variadic_call(std::forward<Args>(args)...);
    }

    //! \brief Enumeration of operations.
    enum class Operation {
        PACKSIZE, //!< Calculating serialization buffer size
        PACK,     //!< Performing serialization
        UNPACK    //!< Performing de-serialization
    };

    //! \brief Predicate for detecting pairs.
    template<class T>
    struct is_pair {
        constexpr static bool value = false;
    };

    template<class T1, class T2>
    struct is_pair<std::pair<T1,T2>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for detecting vectors.
    template<class T>
    struct is_vector {
        constexpr static bool value = false;
    };

    template<class T1, class Allocator>
    struct is_vector<std::vector<T1,Allocator>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for detecting variants.
    template<class T>
    struct is_variant {
        constexpr static bool value = false;
    };

    template<class... Ts>
    struct is_variant<std::variant<Ts...>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for smart pointers.
    template<class T>
    struct is_ptr {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_ptr<std::shared_ptr<T1>> {
        constexpr static bool value = true;
    };

    template<class T1, class Deleter>
    struct is_ptr<std::unique_ptr<T1, Deleter>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for std::optional.
    template<class T>
    struct is_optional {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_optional<std::optional<T1>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for maps
    template<class T>
    struct is_map {
        constexpr static bool value = false;
    };

    template<class Key, class T, class Compare, class Allocator>
    struct is_map<std::map<Key,T,Compare,Allocator>> {
        constexpr static bool value = true;
    };

    template<class Key, class T, class Hash, class KeyEqual, class Allocator>
    struct is_map<std::unordered_map<Key,T,Hash,KeyEqual,Allocator>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for sets
    template<class T>
    struct is_set {
        constexpr static bool value = false;
    };

    template<class Key, class Compare, class Allocator>
    struct is_set<std::set<Key,Compare,Allocator>> {
        constexpr static bool value = true;
    };

    template<class Key, class Hash, class KeyEqual, class Allocator>
    struct is_set<std::unordered_set<Key,Hash,KeyEqual,Allocator>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for arrays
    template<class T>
    struct is_array {
        constexpr static bool value = false;
    };

    template<class T, std::size_t N>
    struct is_array<std::array<T,N>> {
        constexpr static bool value = true;
    };

    //! Detect existence of \c serializeOp member function
    //!
    //! Base case (no \c serializeOp member function)
    template <typename, class = void>
    struct has_serializeOp : public std::false_type {};

    //! Detect existence of \c serializeOp member function
    //!
    //! Non-default, albeit common, case (type has \c serializeOp member
    //! function)
    template <typename T>
    struct has_serializeOp<
        T, std::void_t<decltype(std::declval<T>().serializeOp(std::declval<EclMpiSerializer&>()))>
    > : public std::true_type {};

    //! \brief Handler for pairs.
    template<class T1, class T2>
    void pair(const std::pair<T1,T2>& data)
    {
        if constexpr (has_serializeOp<T1>::value)
            const_cast<T1&>(data.first).serializeOp(*this);
        else
            (*this)(data.first);

        if constexpr (has_serializeOp<T2>::value)
            const_cast<T2&>(data.second).serializeOp(*this);
        else
            (*this)(data.second);
    }

    //! \brief Handler for smart pointers.
    template<class PtrType>
    void ptr(const PtrType& data)
    {
        using T1 = typename PtrType::element_type;
        bool value = data ? true : false;
        (*this)(value);
        if (m_op == Operation::UNPACK && value) {
            const_cast<PtrType&>(data).reset(new T1);
        }
        if (data) {
            if constexpr (has_serializeOp<T1>::value)
                data->serializeOp(*this);
            else
                (*this)(*data);
        }
    }

    Parallel::Communication m_comm; //!< Communicator to broadcast using

    Operation m_op = Operation::PACKSIZE; //!< Current operation
    size_t m_packSize = 0; //!< Required buffer size after PACKSIZE has been done
    int m_position = 0; //!< Current position in buffer
    std::vector<char> m_buffer; //!< Buffer for serialized data
};

}

#endif
