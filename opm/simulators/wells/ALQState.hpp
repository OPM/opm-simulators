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
#include <optional>
#include <string>


namespace Opm {

template<class Scalar>
class ALQState
{
public:
    static ALQState serializationTestObject();

    Scalar get() const;
    void update_default(Scalar value);
    void set(Scalar value);
    bool oscillation() const;
    void update_count(bool increase);
    void reset_count();
    int  get_increment_count() const;
    int  get_decrement_count() const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(current_alq_);
        serializer(default_alq_);
        serializer(alq_increase_count_);
        serializer(alq_decrease_count_);
    }

    bool operator==(const ALQState&) const;

private:
    std::optional<Scalar> current_alq_;
    Scalar default_alq_{0.0};
    int alq_increase_count_{0};
    int alq_decrease_count_{0};
};

}

#endif



