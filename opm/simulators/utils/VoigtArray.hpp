/*
  Copyright 2025 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_UTIL_VOIGT_ARRAY_HPP
#define OPM_UTIL_VOIGT_ARRAY_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <vector>

namespace Opm {

enum class VoigtIndex {
    XX =  0, XY =  5, XZ = 4,
    YX = XY, YY =  1, YZ = 3,
    ZX = XZ, ZY = YZ, ZZ = 2,
};

template<class T>
class VoigtContainer
{
public:
    static constexpr auto indices = std::array{
        Opm::VoigtIndex::XX,
        Opm::VoigtIndex::XY,
        Opm::VoigtIndex::XZ,
        Opm::VoigtIndex::YX,
        Opm::VoigtIndex::YY,
        Opm::VoigtIndex::YZ,
        Opm::VoigtIndex::ZX,
        Opm::VoigtIndex::ZY,
        Opm::VoigtIndex::ZZ,
    };

    static constexpr auto unique_indices = std::array{
        Opm::VoigtIndex::XX,
        Opm::VoigtIndex::YY,
        Opm::VoigtIndex::ZZ,
        Opm::VoigtIndex::YZ,
        Opm::VoigtIndex::XZ,
        Opm::VoigtIndex::XY
    };

    static constexpr auto diag_indices = std::array{
        Opm::VoigtIndex::XX,
        Opm::VoigtIndex::YY,
        Opm::VoigtIndex::ZZ,
    };

    VoigtContainer() = default;

    template<class Array>
    VoigtContainer(const Array& array);

    VoigtContainer(std::initializer_list<T> value)
    {
        std::copy_n(value.begin(),
                    std::min(data_.size(), value.size()),
                    data_.begin());
    }

    const T& operator[](const VoigtIndex idx) const
    { return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)]; }

    T& operator [](const VoigtIndex idx)
    { return data_[static_cast<std::underlying_type_t<VoigtIndex>>(idx)]; }

    constexpr std::size_t size() const { return data_.size(); }

protected:
    std::array<T, 6> data_{};
};

template<class Scalar>
class VoigtArray : public VoigtContainer<std::vector<Scalar>>
{
public:
    VoigtArray() = default;
    explicit VoigtArray(const std::size_t size);

    void resize(const std::size_t size);

    Scalar operator()(const VoigtIndex idx, const std::size_t i) const;
    Scalar& operator()(const VoigtIndex idx, const std::size_t i);

    void assign(const std::size_t i, const VoigtContainer<Scalar>& array);
};

} // namespace Opm

#endif // OPM_UTIL_VOIGT_ARRAY_HPP
