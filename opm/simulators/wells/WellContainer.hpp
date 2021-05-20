/*
  Copyright 2021 Equinor ASA

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

#ifndef OPM_WELL_CONTAINER_HEADER_INCLUDED
#define OPM_WELL_CONTAINER_HEADER_INCLUDED

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {


/*
  The WellContainer<T> class is a small utility class designed to manage the
  dynamic state of per well quantities, like active control and phase rates. The
  values are stored continously in a vector, but they are added with a name, and
  can also be accessed and updated with the name.

  The class is created to facilitate safe and piecewise refactoring of the
  WellStateFullyImplicitBlackOil class, and might have a short life in the
  development timeline.
*/


template <class T>
class WellContainer {
public:

    WellContainer() = default;


    WellContainer(std::initializer_list<std::pair<std::string,T>> init_list) {
        for (const auto& [name, value] : init_list)
            this->add(name, value);
    }


    std::size_t size() const {
        return this->m_data.size();
    }

    void add(const std::string& name, T&& value) {
        if (index_map.count(name) != 0)
            throw std::logic_error("An object with name: " + name + " already exists in container");

        this->index_map.emplace(name, this->m_data.size());
        this->m_data.push_back(std::forward<T>(value));
    }

    void add(const std::string& name, const T& value) {
        if (index_map.count(name) != 0)
            throw std::logic_error("An object with name: " + name + " already exists in container");

        this->index_map.emplace(name, this->m_data.size());
        this->m_data.push_back(value);
    }

    bool has(const std::string& name) const {
        return (index_map.count(name) != 0);
    }


    void update(const std::string& name, T&& value) {
        auto index = this->index_map.at(name);
        this->m_data[index] = std::forward<T>(value);
    }

    void update(const std::string& name, const T& value) {
        auto index = this->index_map.at(name);
        this->m_data[index] = value;
    }

    void update(std::size_t index, T&& value) {
        this->m_data.at(index) = std::forward<T>(value);
    }

    void update(std::size_t index, const T& value) {
        this->m_data.at(index) = value;
    }

    /*
      Will copy the value from other to this - for all wells which are present
      in both containers.
    */
    void copy_welldata(const WellContainer<T>& other) {
        if (this->index_map == other.index_map)
            this->m_data = other.m_data;
        else {
            for (const auto& [name, index] : this->index_map)
                this->update_if(index, name, other);
        }
    }

    /*
      Will copy the value for well @name from other to this. The well @name must
      exist in both containers, otherwise an exception is thrown.
    */
    void copy_welldata(const WellContainer<T>& other, const std::string& name) {
        auto this_index = this->index_map.at(name);
        auto other_index = other.index_map.at(name);
        this->m_data[this_index] = other.m_data[other_index];
    }

    T& operator[](std::size_t index) {
        return this->m_data.at(index);
    }

    const T& operator[](std::size_t index) const {
        return this->m_data.at(index);
    }

    T& operator[](const std::string& name) {
        auto index = this->index_map.at(name);
        return this->m_data[index];
    }

    const T& operator[](const std::string& name) const {
        auto index = this->index_map.at(name);
        return this->m_data[index];
    }

    void clear() {
        this->m_data.clear();
        this->index_map.clear();
    }

    typename std::vector<T>::const_iterator begin() const {
        return this->m_data.begin();
    }

    typename std::vector<T>::const_iterator end() const {
        return this->m_data.end();
    }

    const std::vector<T>& data() const {
        return this->m_data;
    }


private:
    void update_if(std::size_t index, const std::string& name, const WellContainer<T>& other) {
        auto other_iter = other.index_map.find(name);
        if (other_iter == other.index_map.end())
            return;

        auto other_index = other_iter->second;
        this->m_data[index] = other.m_data[other_index];
    }


    std::vector<T> m_data;
    std::unordered_map<std::string, std::size_t> index_map;
};


}


#endif
