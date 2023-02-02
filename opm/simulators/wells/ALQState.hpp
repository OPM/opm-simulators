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

#ifndef OPM_ALQ_STATE_HEADER_INCLUDED
#define OPM_ALQ_STATE_HEADER_INCLUDED

#include <map>
#include <string>
#include <vector>


namespace Opm {

class ALQState {
public:
    static ALQState serializationTestObject();

    std::size_t pack_size() const;
    std::size_t unpack_data(const double * data);
    std::size_t pack_data(double * data) const;

    double get(const std::string& wname) const;
    void update_default(const std::string& wname, double value);
    void set(const std::string& wname, double value);
    bool oscillation(const std::string& wname) const;
    void update_count(const std::string& wname, bool increase);
    void reset_count();
    int  get_increment_count(const std::string& wname) const;
    int  get_decrement_count(const std::string& wname) const;
    void set_debug_counter(int value);
    int  get_debug_counter();
    int  update_debug_counter();

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(current_alq_);
        serializer(default_alq_);
        serializer(alq_increase_count_);
        serializer(alq_decrease_count_);
        serializer(debug_counter_);
    }

    bool operator==(const ALQState&) const;

private:
    std::map<std::string, double> current_alq_;
    std::map<std::string, double> default_alq_;
    std::map<std::string, int> alq_increase_count_;
    std::map<std::string, int> alq_decrease_count_;
    int debug_counter_ = 0;
};


}

#endif



