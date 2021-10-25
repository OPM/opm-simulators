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

#include <dune/common/version.hh>
#include <opm/simulators/utils/ParallelRestart.hpp>

#include <optional>
#include <type_traits>
#include <utility>
#include <variant>

namespace Opm {

/*! \brief Class for (de-)serializing and broadcasting data in parallel.
 *!  \details Can be called on any class with a serializeOp member. Such classes
 *!           are referred to as 'complex types' in the documentation.
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
    //! \tparam complexType Whether or not T is a complex type
    //! \param data The vector to (de-)serialize
    template <typename T, bool complexType = true>
    void vector(std::vector<T>& data)
    {
        auto handle = [&](auto& d)
        {
            for (auto& it : d) {
              if constexpr (is_pair<T>::value)
                  pair(it);
              else if constexpr (is_ptr<T>::value)
                  ptr(it);
              else if constexpr (!complexType)
                  (*this)(it);
              else
                  it.serializeOp(*this);
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

    template <class Array, bool complexType = true>
    void array(Array& data)
    {
        using T = typename Array::value_type;

        auto handle = [&](auto& d) {
            for (auto& it : d) {
                if constexpr (is_pair<T>::value)
                    pair(it);
                else if constexpr (is_ptr<T>::value)
                    ptr(it);
                else if constexpr (!complexType)
                    (*this)(it);
                else
                    it.serializeOp(*this);
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

    //! \brief Handler for std::variant<> with four types

    /*
      The std::variant<> serialization is a first attempt and *not* particularly
      general. In particular that implies:

        1. It is hardcoded to hold exactly four alternative types T0, T1, T2 and
           T3.

        2. All the four types T0, T1, T2 and T3 must implement the ::serializeOp( )
           method. This implies that a variant with a fundamental type like e.g.
           std::variant<int, double, Opm::Well, Opm::Group> will *not* work in the
           current implementation.
    */
    template<class T0, class T1, class T2, class T3>
    void variant(const std::variant<T0,T1,T2,T3>& _data)
    {
        auto handle = [&](auto& d) {
            d.serializeOp(*this);
        };

        std::variant<T0,T1,T2,T3>& data = const_cast<std::variant<T0,T1,T2,T3>&>(_data);
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.index(), m_comm);
            std::visit( [&] (auto& arg) { handle(arg); }, data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.index(), m_buffer, m_position, m_comm);
            std::visit([&](auto& arg) { handle(arg); }, data);
        } else if (m_op == Operation::UNPACK) {
            size_t index;
            Mpi::unpack(index, m_buffer, m_position, m_comm);

            if (index == 0) {
                data = T0();
                handle(std::get<0>(data));
            } else if (index == 1) {
                data = T1();
                handle(std::get<1>(data));
            } else if (index == 2) {
                data = T2();
                handle(std::get<2>(data));
            } else if (index == 3) {
                data = T3();
                handle(std::get<3>(data));
            } else
                std::logic_error("Internal meltdown in std::variant<T0,T1,T2,T3> unpack");
        }
    }


    //! \brief Handler for std::variant<> with two fundamental types

    /*
      This std::variant serialization is highly specialized:

        1. It is hardcoded to take exactly two types T0 and T1.

        2. Both T0 and T1 must be basic types where Mpi::pack(T, ...) overloads
           must exist.

    */
    template<class T0, class T1>
    void variant(const std::variant<T0,T1>& data)
    {
        auto pack_size = [&](auto& d) {
                             m_packSize += Mpi::packSize(d, m_comm);
                         };

        auto packer = [&](auto& d) {
                          Mpi::pack(d, m_buffer, m_position, m_comm);
                      };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.index(), m_comm);
            std::visit( [&] (auto& arg) { pack_size(arg); }, data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.index(), m_buffer, m_position, m_comm);
            std::visit([&](auto& arg) { packer(arg); }, data);
        } else if (m_op == Operation::UNPACK) {
            size_t index;
            std::variant<T0,T1>& mutable_data = const_cast<std::variant<T0,T1>&>(data);
            Mpi::unpack(index, m_buffer, m_position, m_comm);

            if (index == 0) {
                T0 t0;
                Mpi::unpack(t0, m_buffer, m_position, m_comm);
                mutable_data = t0;
            } else if (index == 1) {
                T1 t1;
                Mpi::unpack(t1, m_buffer, m_position, m_comm);
                mutable_data = t1;
            } else
                throw std::logic_error("Internal meltdown in std::variant<T0,T1> unpack loaded index=" + std::to_string(index) + " allowed range: [0,1]");
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
                if constexpr (has_serializeOp<T>::value) {
                    const_cast<T&>(*data).serializeOp(*this);
                } else
                    m_packSize += Mpi::packSize(*data, m_comm);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.has_value(), m_buffer, m_position, m_comm);
            if (data.has_value()) {
                if constexpr (has_serializeOp<T>::value) {
                    const_cast<T&>(*data).serializeOp(*this);
                } else {
                    Mpi::pack(*data, m_buffer, m_position, m_comm);
                }
            }
        } else if (m_op == Operation::UNPACK) {
            bool has;
            Mpi::unpack(has, m_buffer, m_position, m_comm);
            if (has) {
                T res;
                if constexpr (has_serializeOp<T>::value) {
                    res.serializeOp(*this);
                } else {
                    Mpi::unpack(res, m_buffer, m_position, m_comm);
                }
                const_cast<std::optional<T>&>(data) = res;
            }
        }

    }

    //! \brief Handler for maps.
    //! \tparam Map map type
    //! \tparam complexType Whether or not Data in map is a complex type
    //! \param map The map to (de-)serialize
    template<class Map, bool complexType = true>
    void map(Map& data)
    {
        using Key = typename Map::key_type;
        using Data = typename Map::mapped_type;

        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                this->template vector<typename Data::value_type,complexType>(d);
            else if constexpr (is_ptr<Data>::value)
                ptr(d);
            else if constexpr (complexType)
                d.serializeOp(*this);
            else
                (*this)(d);
        };

        auto keyHandle = [&](auto& d)
        {
              if constexpr (is_pair<Key>::value)
                  pair(d);
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

    template<class Set, bool complexType = true>
    void set(Set& data)
    {
        using Data = typename Set::value_type;

        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                this->template vector<typename Data::value_type,complexType>(d);
            else if constexpr (is_ptr<Data>::value)
                ptr(d);
            else if constexpr (complexType)
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
    template<class T>
    void broadcast(T& data)
    {
        if (m_comm.size() == 1)
            return;

        if (m_comm.rank() == 0) {
            try {
                pack(data);
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
            unpack(data);
        }
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

    template<class T1>
    struct is_vector<std::vector<T1>> {
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
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<class T1, class T2>
    void pair(const std::pair<T1,T2>& data)
    {
        if constexpr (std::is_pod<T1>::value || std::is_same<T1,std::string>::value)
            (*this)(data.first);
        else
            data.first.serializeOp(*this);

        if constexpr (std::is_pod<T2>::value || std::is_same<T2,std::string>::value)
            (*this)(data.second);
        else
            const_cast<T2&>(data.second).serializeOp(*this);
    }

    //! \brief Handler for smart pointers.
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<class PtrType>
    void ptr(const PtrType& data)
    {
        using T1 = typename PtrType::element_type;
        bool value = data ? true : false;
        (*this)(value);
        if (m_op == Operation::UNPACK && value) {
            const_cast<PtrType&>(data).reset(new T1);
        }
        if (data)
            data->serializeOp(*this);
    }

    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> m_comm; //!< Communicator to broadcast using

    Operation m_op = Operation::PACKSIZE; //!< Current operation
    size_t m_packSize = 0; //!< Required buffer size after PACKSIZE has been done
    int m_position = 0; //!< Current position in buffer
    std::vector<char> m_buffer; //!< Buffer for serialized data
};

}

#endif
